#include "TrianglePhong.h"
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <iostream>


using namespace std;
using namespace Eigen;

static TrianglePhong::Matrix2 polar2d_fast(const TrianglePhong::Matrix2 &m)
{
  TrianglePhong::Matrix2 ret = TrianglePhong::Matrix2::Zero(2, 2);
  ret(0,0) =  m(1,1);
  ret(1,0) = -m(0,1);
  ret(0,1) = -m(1,0);
  ret(1,1) =  m(0,0);
  if (m.determinant() < 0)
    ret = ret * -1;
  
  ret = ret + m;
  double c0 = sqrt(pow(ret(0,0),2) + pow(ret(1,0),2));
  double c1 = sqrt(pow(ret(0,1),2) + pow(ret(1,1),2));

  ret(0,0) /= c0;
  ret(1,0) /= c0;

  ret(0,1) /= c1;
  ret(1,1) /= c1;

  return ret;
}

static TrianglePhong::Matrix2 polar2d(const TrianglePhong::Matrix2 &m)
{
  JacobiSVD<TrianglePhong::Matrix2> svd(m, ComputeFullU | ComputeFullV);
  TrianglePhong::Matrix2 diag = TrianglePhong::Matrix2::Identity();
  if(m.determinant() < 0)
    diag(1, 1) = -1;
      
  return svd.matrixU() * diag * svd.matrixV().transpose();
}



void TrianglePhong::precompute(Vector corners[3], Basis* tangents, Basis* ert, Basis et1[3])
{
  
  //r's
  Matrix2 r[3];
  for(int i = 0; i < 3; ++i)
    r[i] = polar2d_fast(tangents[(i + 1) % 3] * tangents[i].transpose());
  
  //e's
  Matrix2 e[3];
  for(int i = 0; i < 3; ++i)
    e[i] = Matrix2::Identity();
  
  //all six of these should be close together after multiplication by e.
  Basis b[6];
  for(int i = 0; i < 3; ++i)
  {
    int i1 = (i + 1) % 3;
    b[i * 2] = r[i] * tangents[i];
    b[i * 2 + 1] = tangents[i1];
  }
  
  //iteratively solve for the e's
  for(int times = 0; times < N_ITER; ++times)
  {
    for(int i = 0; i < 3; ++i)
    {
      Matrix2 m = Matrix2::Zero();
      for(int j = 0; j < 3; ++j)
      {
        if(j == i)
          continue;
        for(int i1 = 0; i1 < 2; ++i1)
          for(int j1 = 0; j1 < 2; ++j1)
            m += e[j] * b[j * 2 + j1] * b[i * 2 + i1].transpose();
      }
      
      e[i] = polar2d_fast(m);
    }
    
  }
  
  for(int i = 0; i < 3; ++i)
  {
    ert[i] = e[i] * r[i] * tangents[i];
    et1[i] = e[i] * tangents[(i+1)%3];
  }
}

TrianglePhong::TrianglePhong(MatrixXd Corners, MatrixXd Tangents)
{
  // Alloc aligned memory
  _corners  = igl::aligned_matrix(Corners);
  _tangents = igl::aligned_matrix(Tangents);
  
  _ert = igl::aligned_matrix(Tangents);
  _et1 = igl::aligned_matrix(Tangents);
  
  
  MatrixXf Cornersf  = Corners.cast<float>();
  MatrixXf Tangentsf = Tangents.cast<float>();
  
  // Precompute per triangle
#pragma omp parallel for num_threads(8)
  for(int fid = 0; fid<Corners.rows(); ++fid)
  {
    Vector corners[3];
    Basis  tangents[3];
    Basis  ert[3];
    Basis  et1[3];
    
    for(int i = 0; i < 3; ++i)
    {
      corners[i] = Cornersf.block(fid,i*8,1,8).transpose();
      tangents[i].row(0) = Tangentsf.block(fid,i*16,1,8);
      tangents[i].row(1) = Tangentsf.block(fid,i*16+8,1,8);
    }
    
    precompute(corners, tangents, ert, et1);
    
    // Copy in aligned memory
    for(unsigned i=0; i < 3; ++i)
    {
      float* target_ert = getErt(fid, i);
      for(unsigned j = 0; j < 8; ++j)
      {
        target_ert[j]   = ert[i](0,j);
        target_ert[j+8] = ert[i](1,j);
      }
      
      float* target_et1 = getEt1(fid, i);
      for(unsigned j = 0; j < 8; ++j)
      {
        target_et1[j  ] = et1[i](0,j);
        target_et1[j+8] = et1[i](1,j);
      }
    }
  }
  
  // Precompute per triangle first Newton iteration
  // Alloca a lot of stuff
  _blendB  = igl::aligned_matrix(Corners.rows(),16);
  _blendp  = igl::aligned_matrix(Corners.rows(),8);
  _minusbv = igl::aligned_matrix(Corners.rows(),6);
  _blended = igl::aligned_matrix(Corners.rows(),16*3);
  float w[3];
  for (int i=0; i<3; ++i)
    w[i] = 1./3.;
  float eps = 0.0001;
  for(unsigned fid = 0; fid<Corners.rows(); ++fid)
  {
    float* f_blendB = _blendB.ptr()+fid*16;
    float* f_blendp = _blendp.ptr()+fid*8;
    float* f_minusbv = _minusbv.ptr()+fid*6;
    float* f_blended = _blended.ptr()+fid*16*3;
    
    // Precompute First Iteration
    blendf(fid,w,f_blendB);
    blendPosf(fid,w,f_blendp);
    
    for(int j=0; j<3; ++j)
    {
      const float* v = getCorner(fid, j);
      for(int i=0; i<2; ++i)
      {
        float temp = 0;
        
        for(int t=0; t<8; ++t)
        {
          temp += f_blendB[i*8 + t] * v[t];
        }
        f_minusbv[i*3+j] = -temp;
      }
    }
    
    for(int j = 0; j < 3; ++j)
    {
      float coeff1[6];
      float coeff2[6];
      
      float woffs[3];
      for(int i=0; i<3; ++i)
        woffs[i] = w[i];
      
      woffs[j] += eps;
      computeW(woffs, coeff1);
      woffs[j] -= 2*eps;
      computeW(woffs, coeff2);
      
      for(int i=0; i<6; ++i)
        coeff1[i] = (coeff1[i]-coeff2[i])/(2*eps);
      
      blend_aux(fid, coeff1, &(f_blended[j*16]));
    }
    
  }
 
}

void TrianglePhong::computeW(const float* bary, float coeff[6])
{
  float halfProd = 0.5 * bary[0] * bary[1] * bary[2];
  float scale = bary[0] * bary[1] + bary[1] * bary[2] + bary[0] * bary[2];
  if(fabs(scale) < 1e-8) //avoid division by zero
    scale = 1e-8;
  
  for(int i = 0; i < 3; ++i)
  {
    int i1 = (i + 1) % 3;
    coeff[i*2  ] = (bary[i] * bary[i]  * bary[i1] + halfProd) / scale;
    coeff[i*2+1] = (bary[i] * bary[i1] * bary[i1] + halfProd) / scale;
  }
}

void TrianglePhong::blend_aux(unsigned fid, float coeff[6], float* out)
{
  
  for(unsigned j=0; j<16; ++j)
  {
    out[j] = 0;
  }
  
  for(int i = 0; i < 3; ++i)
  {
    float& coeff1 = coeff[i*2  ];
    float& coeff2 = coeff[i*2+1];
    
    float* ert = getErt(fid,i);
    float* et1 = getEt1(fid,i);
    
#ifdef __INTEL_COMPILER
    __assume_aligned(ert, 32);
    __assume_aligned(et1, 32);
#endif
    
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#pragma vector always
#endif
    for(unsigned j=0; j<16; ++j)
    {
      out[j] += coeff1 * ert[j] + coeff2 * et1[j];
    }
  }
}

void TrianglePhong::blendf(unsigned fid, const float* bary, float* blended)
{
  float coeff[6];
  
  computeW(bary,coeff);
  blend_aux(fid,coeff,blended);
}

void TrianglePhong::blendPosf(unsigned fid, float* bary, float* blended)
{
  const float* c0 = getCorner(fid, 0);
  const float* c1 = getCorner(fid, 1);
  const float* c2 = getCorner(fid, 2);
  
#ifdef __INTEL_COMPILER
  __assume_aligned(c0, 32);
  __assume_aligned(c1, 32);
  __assume_aligned(c2, 32);
#endif
  
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#pragma vector always
#endif
  for(unsigned i=0; i<8; ++i)
    blended[i] = bary[0] * c0[i] + bary[1] * c1[i] + bary[2] * c2[i];
}

void TrianglePhong::project(unsigned fid, const float* point, float* w)
{
  w[0] = 1./3.;
  w[1] = 1./3.;
  w[2] = 1./3.;

  Matrix3 m;
  m.row(2) = RowVector3::Ones();
  
  // Temporary stuff
  //Eigen::Matrix<float,2,3> minusbv;
  float minusbv[6];
  float eps = 0.0001;
  float blended[16];
  float blendB[16];
  float blendp[8];
  float pvb[8];
  float err[2];
  
  // First Iteration
  {
    const float* f_blendB = _blendB.ptr()+fid*16;
    const float* f_blendp = _blendp.ptr()+fid*8;
    const float* f_minusbv = _minusbv.ptr()+fid*6;
    const float* f_blended = _blended.ptr()+fid*16*3;
    
    for(int i=0; i<8; ++i)
      pvb[i] = point[i] - f_blendp[i];
            
    for(int i=0; i<2; ++i)
    {
      float temp = 0;
      for(int t=0; t<8; ++t)
        temp += f_blendB[i*8+t] * pvb[t];
      err[i] = temp;
    }
    
    if ((err[0]*err[0] + err[1]*err[1]) < 1e-12)
      return;

    for(int j = 0; j < 3; ++j)
    {
      
      //m.col(j).head<2>() = blended * pvb + minusbv.col(j);
      for(int i=0; i<2; ++i)
      {
        float temp = 0;
        for(int t=0; t<8; ++t)
          temp += f_blended[j*16+i*8+t] * pvb[t];
        m(i,j) = temp + f_minusbv[i*3+j];
      }
      
    }
    
    

//    cerr << "m :" << m << endl;
//    cerr << "mi:" << m.inverse() << endl;
//
//    cerr << "md :" << md << endl;
//    cerr << "mdi:" << md.inverse() << endl;
//    
//    cerr << "v: " << Vector3(err[0], err[1], 0.) << endl;

    Matrix3d md = m.cast<double>();
    Vector3d wd = md.inverse() * Vector3d(err[0], err[1], 0.);
    
    // w = w - wd;
    for (int i=0;i<3;++i)
      w[i] -= wd[i];
    
    if(w[0] < -0.1 || w[1] < -0.1 || w[2] < -0.1) //anything less than that won't land even after newton
      return;
 
  }
  const int maxNewtonSteps = 10;
  int newtonSteps;
  for(newtonSteps = 0; newtonSteps < maxNewtonSteps; ++newtonSteps)
  {
    blendf(fid,w,blendB);
    blendPosf(fid,w,blendp);
    //Vector pvb = point - blendPos(fid,w);
    
    for(int i=0; i<8; ++i)
      pvb[i] = point[i] - blendp[i];
    
    // Vector2 err = blendB * pvb;
    for(int i=0; i<2; ++i)
    {
      float temp = 0;
      for(int t=0; t<8; ++t)
        temp += blendB[i*8+t] * pvb[t];
      err[i] = temp;
    }
    
    if ((err[0]*err[0] + err[1]*err[1]) < 1e-12)
      break;

    // f(b)  = blend(w(b)) * p - v*b
    // df(b) = blend(dw(b)) * (p-v*b) + blend(b)*(-v)
    
    // Gradient computation

    // minusbv = - blendB * V;

    for(int j=0; j<3; ++j)
    {
      const float* v = getCorner(fid, j);
      for(int i=0; i<2; ++i)
      {
        float temp = 0;

        for(int t=0; t<8; ++t)
        {
          temp += blendB[i*8 + t] * v[t];
        }
        minusbv[i*3+j] = -temp;
      }
    }
    
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#pragma vector always
#endif
    for(int j = 0; j < 3; ++j)
    {
      float coeff1[6];
      float coeff2[6];

      float woffs[3];
      for(int i=0; i<3; ++i)
        woffs[i] = w[i];
      
      woffs[j] += eps;
      computeW(woffs, coeff1);
      woffs[j] -= 2*eps;
      computeW(woffs, coeff2);
      
      for(int i=0; i<6; ++i)
        coeff1[i] = (coeff1[i]-coeff2[i])/(2*eps);
      
      blend_aux(fid, coeff1, blended);
            
      //m.col(j).head<2>() = blended * pvb + minusbv.col(j);
      for(int i=0; i<2; ++i)
      {
        float temp = 0;
        for(int t=0; t<8; ++t)
          temp += blended[i*8+t] * pvb[t];
        m(i,j) = temp + minusbv[i*3+j];
      }
      
    }
    
    Matrix3d md = m.cast<double>();
    Vector3d wd = md.inverse() * Vector3d(err[0], err[1], 0.);
    
    // w = w - wd;
    for (int i=0;i<3;++i)
      w[i] -= wd[i];
    
    if(w[0] < -0.1 || w[1] < -0.1 || w[2] < -0.1) //anything less than that won't land even after newton
      return;
    
  }
  
//  if(newtonSteps == maxNewtonSteps) //may not have converged
//  {
//    cerr << "WTF is happening!" << endl;
//    cerr << point[0] << " " << point[1] << endl;
//    cerr << w[0] << " " << w[1] << " " << w[2] << endl;
//    cerr << "---------------" << endl;
//    //Vector2 err = blend(fid,w) * (point - blendPos(fid,w));
//    
//    if (w[0] != w[0])
//    {
//      cerr << "NAN!" << endl;
//    }
//    
//    if((err[0]*err[0] + err[1]*err[1]) > 1e-12) //check the error
//    {
//      //return Bary(-1, -1, 3); //invalid if newton didn't converge
//      w[0] = -1;
//      w[1] = -1;
//      w[2] = 3;
//      return;
//    }
//  }
}
