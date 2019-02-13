#include "Phong.h"

#include "igl/tt.h"
#include "igl/vf.h"
#include "igl/edgetopology.h"
#include "igl/per_vertex_normals.h"
#include "igl/Timer.h"
#include "igl/adjacency_list.h"

#include <float.h>
#include <tr1/unordered_set>
#include "null.h"


// get AVX intrinsics  
#ifdef ICC
#include <immintrin.h>  
static const __m128 SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
static const __m128 eye = _mm_setr_ps(1, 0, 0, 0);
static const __m128 zero = _mm_setzero_ps();
static const __m128 ones = _mm_set1_ps(1);
#endif

//#define TIMER_PHONG
#include <queue>

//#define TIME_BACKWARD

//#define REGISTERS_ONLY

using namespace std;

Anchor_Set_Phong::Anchor_Set_Phong(const std::vector<PointOnSurface > &point_set, const MatrixXX& E, const MatrixX3i& F):
_point_set(point_set)
{
  anchors_E.resize(point_set.size(),8);
  for(unsigned i=0; i<point_set.size();++i)
    anchors_E.row(i) = point_set[i].get_position(E, F);
  
  a_anchors_E = igl::aligned_matrix(anchors_E);
}

Anchor_Set_Phong::~Anchor_Set_Phong()
{
  
}


Phong::Phong()
{
}

Phong::~Phong()
{
  a_V.dealloc();
}

void Phong::init(const MatrixXX& V, const MatrixXXi& F, vector<Basis>* basis, bool x_bUseEuclideanInit)
{
  _V = V;
  _F = F;
  
  // AVX aligned matrices
  a_V = igl::aligned_matrix(_V);
  a_V_cols = _V.cols();
  
  if (basis)
    _basis = *basis;
  else
    if (x_bUseEuclideanInit)
      _basis = computeBasisE(_V,_F);
    else
      _basis = computeBasis(_V,_F);
  
  _a_basis.resize(_basis.size());
  for (size_t i = 0 ; i <_basis.size(); ++i)
    _a_basis[i] = igl::aligned_matrix(_basis[i]);
  
  // Compute topological relations
  igl::adjacency_list(_F, VV);
  igl::tt(_V, _F, TT, TTi);
  igl::vf(V, F, VF, VFi);
  
  precomputeBarycentricOperator();
  precomputeKRings();

  // Precomputation
  Eigen::MatrixXd Corners(_F.rows(),24);
  Eigen::MatrixXd Tangents(_F.rows(),48);
  
  for (unsigned i=0; i<_F.rows(); ++i)
  {
    for (unsigned j=0;j<3;++j)
      Corners.block(i,j*8,1,8) = _V.row(_F(i,j));
    for (unsigned j=0;j<3;++j)
    {
      Tangents.block(i,j*16,1,8)   = _basis[_F(i,j)].row(0);
      Tangents.block(i,j*16+8,1,8) = _basis[_F(i,j)].row(1);
    }
  }

  TrianglesPhong = TrianglePhong(Corners,Tangents);

}

void Phong::precomputeBarycentricOperator()
{
  _baryop.resize(_F.rows(),3);
  _baryop_cropped.resize(_F.rows(),3);
  
  MatrixXX _baryop_col_07(_F.rows(),3*3*8);
  MatrixXX _baryop_col_8(_F.rows(),3*3);
  MatrixXX _baryop_col_8_f(_F.rows(),3*4);
  _baryop_col_8_f.setZero();
  
  for(unsigned fid=0; fid<_F.rows();++fid)
  {
    // Positions of the triangle vertices
    Vector8 v[3]; 
    for (unsigned i=0; i<3; ++i)
      v[i] = _V.row(_F(fid,i));
    
    // Basis of the tangent plane at every vertex
    Basis b[3];    
    for (unsigned i=0; i<3; ++i)
      b[i] = _basis[_F(fid,i)];
    
    for(int i = 0; i < 3; ++i)
    {
      // Bary.inverse() is the projection operator
      Matrix33 bary;
      bary.row(2).setOnes();
      for(int j = 0; j < 3; ++j)
        bary.col(j).head<2>() = b[i] * v[j];
      
      // Augmented basis
      Eigen::Matrix<ScalarType,3,9> augbasis;
      augbasis.setZero(); 
      augbasis.block<2,8>(0,0) = b[i];
      augbasis(2,8) = 1;
      
      _baryop(fid,i) = bary.inverse()*augbasis;
      _baryop_cropped(fid,i) = _baryop(fid,i).block<2,9>(0, 0);
      
      _baryop_col_07.block<1,3*8>(fid,i*8*3) << 
      _baryop(fid,i).block<1,8>(0, 0),
      _baryop(fid,i).block<1,8>(1, 0),
      _baryop(fid,i).block<1,8>(2, 0);
      
      _baryop_col_8.block<1,3>(fid,i*3) << 
      _baryop(fid,i)(0, 8),
      _baryop(fid,i)(1, 8),
      _baryop(fid,i)(2, 8);
      
      _baryop_col_8_f(fid,0*4+i) = _baryop(fid,i)(0, 8);
      _baryop_col_8_f(fid,1*4+i) = _baryop(fid,i)(1, 8);
      _baryop_col_8_f(fid,2*4+i) = _baryop(fid,i)(2, 8);
      
    }
  }
  
  _a_baryop_col_07 = igl::aligned_matrix(_baryop_col_07);
  _a_baryop_col_8 = igl::aligned_matrix(_baryop_col_8);
  _a_baryop_col_8_f = igl::aligned_matrix(_baryop_col_8_f);
}


vector<Phong::Basis> Phong::computeBasis(const MatrixXX& V, const MatrixXXi& F)
{
  vector<vector<int> > VV;
  igl::adjacency_list(F, VV,true);
  
  vector<Basis> basis(V.rows());
  
  for (unsigned i=0; i<V.rows(); ++i)
  {
    RowVectorX t1 = RowVectorX::Zero(V.cols());
    RowVectorX t2 = RowVectorX::Zero(V.cols());
    
    int N = VV[i].size();
    
    for(unsigned j=0; j<N; ++j)
    {
      double w1 = cos((2.0*M_PI*(double)j)/((double)N));
      double w2 = sin((2.0*M_PI*(double)j)/((double)N));
      
      t1 += V.row(VV[i][j]) * w1;
      t2 += V.row(VV[i][j]) * w2;
    }
    
    t1.normalize();
    t2 = t2 - t2.dot(t1)*t1;
    t2.normalize();
    
    Basis b(2,V.cols());
    b.row(0) = t1;
    b.row(1) = t2;
    
    basis[i] = b;      
    
  }
  
  return basis;
}

vector<Phong::Basis> Phong::computeBasisE(const MatrixXX& V, const MatrixXXi& F)
{
  vector<Basis> basis(V.rows());
  MatrixXX V3 = V.block(0, 0, V.rows(), 3);

  MatrixXX N;
  igl::per_vertex_normals(V3, F, N);
  
  for (unsigned i=0; i<V.rows(); ++i)
  {
    
    MatrixXX TP;
    MatrixXX r = N.row(i);
    nullspace(r,TP);
    TP.transposeInPlace();

    RowVectorX t1 = TP.row(0);
    RowVectorX t2 = TP.row(1);
    
    t1.normalize();
    t2 = t2 - t2.dot(t1)*t1;
    t2.normalize();
    
    Basis b = Basis::Zero(2,V.cols());
    b.row(0).head(3) = t1;
    b.row(1).head(3) = t2;
    
    basis[i] = b;      
    
  }
  
  return basis;
}

bool Phong::projectOnTriangle(const unsigned &fid, const Vector8 &p, RowVector3& w)
{  
  // Every row stores the barycentric coordinates of the projections wrt to a tangent plane
  Matrix33 C;
  
  for(int i = 0; i < 3; ++i)
  {
    C.row(i) = _baryop(fid,i).block<3,8>(0,0) * p;
    C(i,0) += _baryop(fid,i)(0,8);
    C(i,1) += _baryop(fid,i)(1,8);
    C(i,2) += _baryop(fid,i)(2,8);
    
  }
  
  // Find w that satisfies C'w = 1*w (eigenvector of w with associated ev 1)
  C -= Matrix33::Identity();
  
  Vector3 w1 = C.col(0).cross(C.col(1));
  Vector3 w2 = C.col(0).cross(C.col(2));
  Vector3 w3 = C.col(1).cross(C.col(2));
  
  w = w1;
  if(fabs(w2.sum()) > fabs(w.sum()))
    w = w2;
  if(fabs(w3.sum()) > fabs(w.sum()))
    w = w3;
  
  w /= w.sum();
  
  return (w(0) > 0) && (w(1) > 0) && (w(2) > 0); 
}
#ifdef ICC

#ifdef REGISTERS_ONLY
//1 cross product at a time (still SSE-d), as much in registers as possible
bool Phong::projectOnTriangle_fast(const unsigned &fid,
                                   const float *p,
                                   float *w)
{
  __m256 avx_p = _mm256_load_ps(p); 
  __m256 avx_ptr, avx_t;
  
  float *ptr =_a_baryop_col_07.ptr()+3*3*8*fid;
  float *ptr8 = _a_baryop_col_8_f.ptr() + 3*4*fid;
  
  __m128 avx_c0, avx_c1, avx_c2;//contain columns 0, 1, 2 of C
  
  //  construct matrix C
  //row 0
  //element 00
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 0 0 x 0 0 0 y] where C00 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF1);
  // avx_c0_: [0 0 0 C00] 
  avx_c0 = _mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1));
  avx_c0 = _mm_sub_ps(avx_c0,eye);
  //element 01
  ptr += 8;
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 0 0 x 0 0 0 y] where C01 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF1);
  // avx_c1_: [0 0 0 C01] 
  avx_c1 = _mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1));
  //element 02
  ptr += 8;
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 0 0 x 0 0 0 y] where C02 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF1);
  // avx_c2_: [0 0 0 C02] 
  avx_c2 = _mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1));
  
  //row 1
  //element 10
  ptr += 8;
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 0 x 0 0 0 y 0] where C10 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF2);
  // avx_c0_: [0 0 C10 C00] 
  avx_c0 = _mm_add_ps(avx_c0,_mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1)));
  //element 11
  ptr += 8;
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 0 x 0 0 0 y 0] where C11 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF2);
  // avx_c1_: [0 0 C11 C01] 
  avx_c1 = _mm_add_ps(avx_c1,_mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1)));
  avx_c1 = _mm_sub_ps(avx_c1,_mm_shuffle_ps(eye, eye, _MM_SHUFFLE(3, 2, 0, 1)));
  //element 12
  ptr += 8;
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 0 x 0 0 0 y 0] where C12 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF2);
  // avx_c2_: [0 0 C12 C02] 
  avx_c2 = _mm_add_ps(avx_c2,_mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1)));
  
  //row 2
  //element 20
  ptr += 8;
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 x 0 0 0 y 0 0] where C20 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF4);
  // avx_c0_: [0 C20 C10 C00] 
  avx_c0 = _mm_add_ps(avx_c0,_mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1)));
  //element 21
  ptr += 8;
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 x 0 0 0 y 0 0] where C21 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF4);
  // avx_c1_: [0 C21 C11 C01] 
  avx_c1 = _mm_add_ps(avx_c1,_mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1)));
  //element 22
  ptr += 8;
  avx_ptr = _mm256_load_ps(ptr); 
  // avx_t: [0 x 0 0 0 y 0 0] where C22 = x+y 
  avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF4);
  // avx_c2_: [0 C22 C12 C02] 
  avx_c2 = _mm_add_ps(avx_c2,_mm_add_ps(_mm256_extractf128_ps(avx_t,0), _mm256_extractf128_ps(avx_t,1)));
  avx_c2 = _mm_sub_ps(avx_c2,_mm_shuffle_ps(eye, eye, _MM_SHUFFLE(3, 0, 1, 2)));
  
  avx_c0 = _mm_add_ps(avx_c0,_mm_load_ps(ptr8) );
  avx_c1 = _mm_add_ps(avx_c1,_mm_load_ps(ptr8+4) );
  avx_c2 = _mm_add_ps(avx_c2,_mm_load_ps(ptr8+8) );
  
  __m128 avx_w1 = CrossProduct(avx_c0,avx_c1);
  __m128 avx_w2 = CrossProduct(avx_c0,avx_c2);
  __m128 avx_w3 = CrossProduct(avx_c1,avx_c2);
  
  avx_c1 = _mm_hadd_ps (zero, avx_w3);//0, 0, w3[2], w3[1]+w3[0] 
  avx_c0 = _mm_hadd_ps (avx_w1, avx_w2);//w1[0]+w1[1], w1[2]+0, w2[0]+w2[1], w2[2]+0 
  avx_c1 = _mm_hadd_ps (avx_c0, avx_c1);//w1[0]+w1[1]+w1[2]+0, w2[0]+w2[1]+w2[2]+0, 0, w3[0]+w3[1]+w3[2]
  
  avx_c2 = _mm_andnot_ps(SIGNMASK, avx_c1);
  
  //compare [fs1 fs2 0 fs3] with [fs2 fs3 0 fs1]
  int comp_res = _mm_movemask_ps(_mm_cmpgt_ps (avx_c2, _mm_shuffle_ps(avx_c2, avx_c2, _MM_SHUFFLE(2, 0, 1, 3)) ) ); 
  switch (comp_res) 
  {
    case 12:
    case 8:
      //fs1 is largest
      avx_c0 = _mm_mul_ps (avx_w1, _mm_div_ps(ones,_mm_shuffle_ps(avx_c1, avx_c1, _MM_SHUFFLE(3, 3, 3, 3))));
      break;
    case 4:
    case 5:
      //fs2 is largest
      avx_c0 = _mm_mul_ps (avx_w2, _mm_div_ps(ones,_mm_shuffle_ps(avx_c1, avx_c1, _MM_SHUFFLE(2, 2, 2, 2))));
      break;
    case 1:
    case 9:
      //fs3 is largest
      avx_c0 = _mm_mul_ps (avx_w3, _mm_div_ps(ones,_mm_shuffle_ps(avx_c1, avx_c1, _MM_SHUFFLE(0, 0, 0, 0))));
      break;
      
    default:
      cerr<<"no match"<<endl;
      break;
  }
  
  _mm_store_ps (w,avx_c0);
  
  return _mm_movemask_ps(_mm_and_ps(avx_c0,SIGNMASK));
}
#else
//1 cross product at a time (still SSE-d), memory and registers
bool Phong::projectOnTriangle_fast(const unsigned &fid,
                                   const float *p,
                                   float *w)
{
  if (!gapfiller)
  {
    __m256 avx_p = _mm256_load_ps(p);
    __m256 avx_ptr, avx_t;
    float f[8];
    float C[16]; // col major
    for(int i = 0; i < 3; ++i)
    {
      const float *ptr =  _a_baryop_col_07.ptr()+3*3*8*fid+ 3*8*i;
      for(int col = 0; col < 3; ++col)
      {
        avx_ptr = _mm256_load_ps(ptr+8*col);
        avx_t = _mm256_dp_ps(avx_ptr,avx_p,0xF1);
        _mm256_store_ps(f,avx_t);
        C[col*4+i] = (f[0]+f[4]);
        C[col*4+i] += _a_baryop_col_8.ptr()[3*3*fid + 3*i + col];
      }
      C[i*4+3] =0.;
      
      C[i*4+i] -=1.;
      
    }
    for(int i = 0; i < 4; ++i)
      C[3*4+i] =0.;
    
    float w1[4], w2[4], w3[4];
    __m128 avx_c0 = _mm_setr_ps(C[0], C[1], C[2], 0);
    __m128 avx_c1 = _mm_setr_ps(C[4], C[5], C[6], 0);
    __m128 avx_c2 = _mm_setr_ps(C[8], C[9], C[10], 0);
    __m128 avx_w1 = CrossProduct(avx_c0,avx_c1);
    __m128 avx_w2 = CrossProduct(avx_c0,avx_c2);
    __m128 avx_w3 = CrossProduct(avx_c1,avx_c2);
    
    _mm_store_ps (w1, avx_w1);
    _mm_store_ps (w2, avx_w2);
    _mm_store_ps (w3, avx_w3);
    
    avx_c0 = _mm_hadd_ps (avx_w1, avx_w2);//w1[0]+w1[1], w1[2]+0, w2[0]+w2[1], w2[2]+0
    avx_c0 = _mm_hadd_ps (avx_c0, avx_w3);//w1[0]+w1[1]+w1[2]+0, w2[0]+w2[1]+w2[2]+0, w3[0]+w3[1], w3[2]+0
    _mm_store_ps (f, avx_c0);
    float s1 = f[0];
    float s2 = f[1];
    float s3 = f[2]+f[3];
    
    float fs1 = fabs(s1);
    float fs2 = fabs(s2);
    float fs3 = fabs(s3);
    
    if(fs1>fs2)
    {
      if(fs1>fs3)
      {
        s1 = 1./s1;
        w[0] = w1[0]*s1;
        w[1] = w1[1]*s1;
        w[2] = w1[2]*s1;
      }
      else
      {
        s3 = 1./s3;
        w[0] = w3[0]*s3;
        w[1] = w3[1]*s3;
        w[2] = w3[2]*s3;
      }
    }
    else
    {
      if(fs2>fs3)
      {
        s2 = 1./s2;
        w[0] = w2[0]*s2;
        w[1] = w2[1]*s2;
        w[2] = w2[2]*s2;
      }
      else
      {
        s3 = 1./s3;
        w[0] = w3[0]*s3;
        w[1] = w3[1]*s3;
        w[2] = w3[2]*s3;
      }
      
    }
  }
  else
  {
    TrianglesPhong.project(fid,p,w);
  }

return (w[0] > EPSNEG) && (w[1] > EPSNEG) && (w[2] > EPSNEG);
}
#endif

#else
//non-AVX

bool Phong::projectOnTriangle_fast(const unsigned &fid,
                                   const float *p,
                                   float *w)
{
  TrianglesPhong.project(fid,p,w);

  return (w[0] >= EPSNEG) && (w[1] >= EPSNEG) && (w[2] >= EPSNEG);  
}

#endif

bool Phong::projectOnTriangleEuclidean_fast(const unsigned &fid,
                                   const float *p,
                                   float *w)
{
  // Positions of the triangle vertices
  Vector8 v[3]; 
  for (unsigned i=0; i<3; ++i)
    v[i] = _V.row(_F(fid,i));
    
  // Basis of the tangent plane at the triangle
  
  Vector8 e0 = (v[1] - v[0]).normalized();
  Vector8 e1 = (v[2] - v[0]).normalized();
  
  Basis b;    

  b.row(0) = e0;
  double s = e0.transpose() * e1;
  b.row(1) = (e1 - s * e0).normalized();
  
  // Bary.inverse() is the projection operator
  Matrix33 bary;
  bary.row(2).setOnes();
  for(int j = 0; j < 3; ++j)
    bary.col(j).head<2>() = b * v[j];

  Vector8 pv;
  for (unsigned i=0;i<8;++i)
    pv(i) = p[i];
  
  Vector3 pp;
  pp.head(2) = b * pv;
  pp(2) = 1;
  
  Vector3 wv; 
  wv = bary.inverse() * pp;
  
  w[0] = wv(0);
  w[1] = wv(1);
  w[2] = wv(2);
  
  return (w[0] > -0.001) && (w[1] > -0.001) && (w[2] > -0.001); 
}

bool Phong::projectBruteForce(const Vector8 &p, int& fid, RowVector3& w, bool Phong)
{
  fid = -1;
  float dist = FLT_MAX;

  float p_[8];
  p_[0]= p[0];
  p_[1]= p[1];
  p_[2]= p[2];
  p_[3]= p[3];
  p_[4]= p[4];
  p_[5]= p[5];
  p_[6]= p[6];
  p_[7]= p[7];
  float w_[3];
  Vector3 cw;

  // Project inside triangles
  for(unsigned i=0;i<_F.rows();++i)
  {
    bool b;
    
    if (Phong)
      b = projectOnTriangle_fast(i, p_, w_);
    else
      b = projectOnTriangleEuclidean_fast(i, p_, w_);

    if (b)
    {
      cw[0] = w_[0];
      cw[1] = w_[1];
      cw[2] = w_[2];
      PointOnSurface ps(i, cw);
      Vector8 pos = ps.get_position(_V,_F);
      
      double currdist = (pos-p).norm();
      
      if (currdist < dist)
      {
        fid = i;
        dist = currdist;
        w = cw;
      }
    }
  }

  if (!Phong)
  {
    // Find closest vertex
    for(unsigned i=0;i<_F.rows();++i)
    {
      for(unsigned j=0;j<3;++j)
      {
        double currdist = (_V.row(_F(i,j)).transpose()-p).norm();
        if (currdist < dist)
        {
          fid = i;
          dist = currdist;
          
          w = RowVector3::Zero();
          w(j) = 1;
        }
      }
    }
    
    // Find closest point on edges
    for(unsigned i=0;i<_F.rows();++i)
    {
      for(unsigned j=0;j<3;++j)
      {
        Vector8 v1 = _V.row(_F(i,j));
        Vector8 v2 = _V.row(_F(i,(j+1)%3));
        
        Vector8 e = v2-v1;
        Vector8 dir = e.normalized();
        
        double t = (p-v1).transpose()*dir;
        t = t / e.norm();
        
        if (t>=0 && t<=1)
        {
          Vector8 pp =  v1 + (t*e.norm()) * dir;
          
          double currdist = (pp-p).norm();
          
          if (currdist < dist)
          {
            fid = i;
            dist = currdist;
            
            w = RowVector3::Zero();
            w(j) = 1-t;
            w((j+1)%3) = t;
          }
        }
      }
    }
    
    
  }
  
  return fid!=-1;
}

bool Phong::project(const Vector8 &p, int fid_start, int& fid, RowVector3& w)
{
  
  // Find closest point
  //WARM START IS A FACE
  int t_vid = _F(fid_start,0);
  t_vid = findClosest(p, t_vid);
  fid_start = VF[t_vid][0];
  
  //  fid_start = findClosestFace(p, fid_start);
  
  float p_[8];
  p_[0]= p[0];
  p_[1]= p[1];
  p_[2]= p[2];
  p_[3]= p[3];
  p_[4]= p[4];
  p_[5]= p[5];
  p_[6]= p[6];
  p_[7]= p[7];
  float w_[3];
  
  
  // Project on the 1-ring
  int t_fid;
  for(unsigned i = 0; i<KRings[fid_start].size(); ++i)
  {
    t_fid = KRings[fid_start][i];
    
    if (projectOnTriangle_fast(t_fid, p_, w_))
    {
      fid = t_fid;
      w[0] = w_[0];
      w[1] = w_[1];
      w[2] = w_[2];
      return true;
    }
  }
  
  //t_fid or fid_start here?
  if (projectBFS(p_, t_fid, fid, w_)) 
  {
    w[0] = w_[0];
    w[1] = w_[1];
    w[2] = w_[2];
    return true;
  }
  else
    cerr<<"projectBFS did not finish"<<endl;
  return false;
}

void Phong::unproject_fast(const Anchor_Set_Phong& anchor_set,
                           const PointOnSurface &point,
                           RowVectorX &weights,
                           const double &ilya_hack_parameter)
{  
#ifdef TIME_BACKWARD
  
  igl::Timer timer;
  double t_init=0, t_sqd=0, t_C = 0, t_solve = 0;
  for (int iteration = 0; iteration <1000; ++iteration)
  {
    
    
#endif
    
#ifdef TIME_BACKWARD
    timer.start();
#endif
    
    const float *bary = point.bc_ptr;
        
    int num_anchors = anchor_set.anchors_E.rows();
#ifdef TIME_BACKWARD
    timer.stop();
    t_init += timer.getElapsedTimeInMicroSec();
#endif
    
#ifdef TIME_BACKWARD
    timer.start();
#endif

    // Computing square distances
    VectorX sqd(num_anchors);

    float baryT[3];
    for (unsigned i=0; i<3; ++i)
      baryT[i] = bary[i];
    
    for (unsigned i=0; i<3; ++i)
      if (bary[i] == 1)
      {
        baryT[i]       = 0.9998;
        baryT[(i+1)%3] = 0.0001;
        baryT[(i+2)%3] = 0.0001;
      }
    
    float p[8];
    TrianglesPhong.blendPosf(point.fid,baryT,p);
    
    float basis[16];
    TrianglesPhong.blendf(point.fid,baryT,basis);

#define SQDNEW
#ifdef SQDNEW
    // Compute tangent plane
          
      float temp2[2];
      float vT[8];
      float vN[8];
      
      // compute vT and vP
      for(int a=0; a<num_anchors; ++a)
      {
        // v = anchors
        float v[8];
        for (int j=0; j<8; ++j)
        {
          v[j] = anchor_set.a_anchors_E.ptr()[a*8+j] - p[j];
        }
        
        // temp2 = T * v
        for(int i=0; i<2; ++i)
        {
          float t = 0;
          for(int j=0; j<8;++j)
            t += basis[i*8 + j] * v[j];
          temp2[i] = t;
        }
        
        // vT = T' * temp2
        for(int j=0; j<8; ++j)
          vT[j] = basis[0+j] * temp2[0] + basis[8+j] * temp2[1];

        // vN = v - vT
        for(int j=0; j<8; ++j)
          vN[j] = v[j] - vT[j];
        
        sqd(a) = 0;
        
        for(int j=0; j<8; ++j)
          sqd(a) += vT[j] * vT[j] + ilya_hack_parameter * vN[j] * vN[j];
      }
        
#else
    
        if (!(point.has_e))
        {
          igl::aligned_matrix a_v((MatrixXX(3,8) << _V.row(_F(point.fid,0)), _V.row(_F(point.fid,1)), _V.row(_F(point.fid,2)) ).finished());
          PointOnSurface::evaluate_average(bary, 3, a_v.ptr(), tt);
          a_v.dealloc();
    
        }
        const float *a_target = (point.has_e)? point.e: tt;
    
    // Compute the gradients using strong distributivity
    //#ifndef ICC
    float vec[8];
    float tx, ty;
    float temp;
    float t;
    for(int i=0; i<num_anchors; ++i)
    {
      temp = 0;
      t = 0;
      for (int j =0; j<8; ++j)
      {
        temp = anchor_set.a_anchors_E.ptr()[i*8+j] - a_target[j];
        vec[j] = temp;
        t += temp * temp;
      }
      sqd(i) = (1 + ilya_hack_parameter)*t;
      
      temp = 0;
      tx = 0;
      for (int j =0; j<8; ++j)
        tx += _a_basis[_F(point.fid,0)].ptr()[0*8+j]*vec[j];
      ty = 0;
      for (int j =0; j<8; ++j)
        ty += _a_basis[_F(point.fid,0)].ptr()[1*8+j]*vec[j];
      temp += bary[0] * (tx*tx + ty*ty );
      tx = 0;
      for (int j =0; j<8; ++j)
        tx += _a_basis[_F(point.fid,1)].ptr()[0*8+j]*vec[j];
      ty = 0;
      for (int j =0; j<8; ++j)
        ty += _a_basis[_F(point.fid,1)].ptr()[1*8+j]*vec[j];
      temp += bary[1] * (tx*tx + ty*ty );
      tx = 0;
      for (int j =0; j<8; ++j)
        tx += _a_basis[_F(point.fid,2)].ptr()[0*8+j]*vec[j];
      ty = 0;
      for (int j =0; j<8; ++j)
        ty += _a_basis[_F(point.fid,2)].ptr()[1*8+j]*vec[j];
      temp += bary[2] * (tx*tx + ty*ty );
      sqd(i) -= ilya_hack_parameter * temp;
    }

#endif
    
#ifdef TIME_BACKWARD
    timer.stop();
    t_sqd += timer.getElapsedTimeInMicroSec();
#endif
    
#ifdef TIME_BACKWARD
    timer.start();
#endif

    // Build the constraints matrix
    float C_0[8], C_1[8];
#ifndef ICC
    for (int l =0; l<8; ++l)
    {
      C_0[l] =0;
      C_1[l] =0;
      for (int k = 0; k<3; ++k)
      {
        C_0[l] += bary[k] * _a_baryop_col_07.ptr()[3*3*8*point.fid+ 3*8*k + 8*0 + l];
        C_1[l] += bary[k] * _a_baryop_col_07.ptr()[3*3*8*point.fid+ 3*8*k + 8*1 + l];
      }
    }
#else
    __m256 avx_bop0, avx_bop1, avx_w, avx_t0,avx_t1;
    __m256 avx_acc0 = _mm256_setzero_ps();
    __m256 avx_acc1 = _mm256_setzero_ps();
    for (int k = 0; k<3; ++k)
    {
      avx_bop0  = _mm256_load_ps(_a_baryop_col_07.ptr()+3*3*8*point.fid+ 3*8*k + 8*0);
      avx_bop1  = _mm256_load_ps(_a_baryop_col_07.ptr()+3*3*8*point.fid+ 3*8*k + 8*1);
      avx_w  = _mm256_set1_ps(bary[k]);
      
      avx_t0   = _mm256_mul_ps(avx_bop0,avx_w);
      avx_acc0    = _mm256_add_ps(avx_t0,avx_acc0);
      avx_t1   = _mm256_mul_ps(avx_bop1,avx_w);
      avx_acc1    = _mm256_add_ps(avx_t1,avx_acc1);
    }
    _mm256_store_ps(C_0,avx_acc0);
    _mm256_store_ps(C_1,avx_acc1);
    
#endif
    
    MatrixXX C (3, num_anchors);
    Vector3 bary_for_solve;

    {
      float* anchors = anchor_set.a_anchors_E.ptr();
      
      float temp2;
      for(int i=0; i<2; ++i)
      {
        for(int j=0; j<num_anchors; ++j)
        {
          temp2=0;
          for(int t=0;t<8;++t)
            temp2 += basis[i*8+t] * anchors[j*8+t];
          C(i,j) = temp2;
        }
      }
      
      C.row(2).setOnes();
      
#ifdef TIME_BACKWARD
      timer.stop();
      t_C+= timer.getElapsedTimeInMicroSec();
#endif
      
#ifdef TIME_BACKWARD
      timer.start();
#endif
      float pp0 = 0;
      float pp1 = 0;
      for(int i=0; i<8; ++i)
      {
        pp0 += basis[i] * p[i];
        pp1 += basis[8+i] * p[i];
      }
      
      bary_for_solve[0] = pp0;
      bary_for_solve[1] = pp1;
      bary_for_solve[2] = 1.;
    }
    
    // Build the energy matrix
    Eigen::DiagonalMatrix<ScalarType,Eigen::Dynamic> Dinv = sqd.cwiseInverse().asDiagonal();
    MatrixXX CDinv  = C*Dinv;
    Matrix33 CCDD = CDinv * CDinv.transpose();
#ifndef TIME_BACKWARD
    if(fabs(CCDD.determinant()) > 1e-7)
    {
#endif
      weights = Dinv * CDinv.transpose() * CCDD.ldlt().solve(bary_for_solve);
#ifndef TIME_BACKWARD
    }
    else
    {
      Eigen::JacobiSVD<MatrixXX> svd(CDinv, Eigen::ComputeFullU | Eigen::ComputeFullV);
      weights = (Dinv*svd.solve(bary_for_solve)).transpose();
    }
#endif
#ifdef TIME_BACKWARD
    timer.stop();
    t_solve += timer.getElapsedTimeInMicroSec();
#endif
#ifdef TIME_BACKWARD
  }
  
  cerr<<"unproject -- init (us): "<<t_init<<endl;
  cerr<<"unproject -- sqd (us): "<<t_sqd<<endl;
  cerr<<"unproject -- C (us): "<<t_C<<endl;
  cerr<<"unproject -- solve (us): "<<t_solve<<endl;
#endif
}


bool Phong::projectBFS(const float *p, int fid_start, int& fid, float* w)
{
  std::vector<int> Q;
  Q.reserve(25);
  
  std::tr1::unordered_set<int> visited;
  size_t iter = 0;
  
  Q.push_back(fid_start);
  visited.insert(fid_start);
  while(iter<Q.size())
  { 
    fid = Q[iter++];
    
    if (projectOnTriangle_fast(fid, p, w))
      return true;
    
    for (int i = 0; i<3; ++i)
    { 
      int next = TT(fid,i);
      if(next >=0 &&  visited.insert(next).second) 
        Q.push_back(next); 
    }
  } 
  
  return false;
}

bool Phong::findMinimumJumpNeg(const float *p, int fid_start, int& fid, float* w)
{
  
  int current_i = fid_start;
  bool b;
  int index_to_jump;
  //  float least_negative_weight;
  float most_negative_weight;
  std::tr1::unordered_set<int> visited;
  bool foundVisited = false;
  
  while (!foundVisited)
  {
    std::pair<std::tr1::unordered_set<int>::iterator,bool> ir = visited.insert(current_i);
    if (!ir.second)
    {
      foundVisited = true;
      cerr<<"Found already visited face. Returning.."<<endl;
      break;
    }
    b = projectOnTriangle_fast(current_i, p, w);
    if (b)
    {
      fid = current_i;
      return true;
    }
    
    index_to_jump = -1;
    //    least_negative_weight = -1e5;
    most_negative_weight = 1e5;
    for (int i = 0; i<3; ++i)
    {
      //      if (w[i] <0 && w[i]>least_negative_weight)
      if (w[i]<most_negative_weight)
      {
        index_to_jump = TT(current_i,i);
        if (index_to_jump >=0)
          //          least_negative_weight = w[i];
          most_negative_weight = w[i];
      }
    }
    if (index_to_jump >=0)
    {
      current_i = index_to_jump;
    }
  }
  return false;
}

int Phong::findClosest(const Vector8& p, int vid)
{
  igl::aligned_matrix a_p((RowVector8)(p.transpose()));
  
  float current_v = squaredNorm(&(a_V.ptr()[vid*a_V_cols]),a_p.ptr());
  int current_i = vid;
  int previous_i = -1;
  
  while (current_i != previous_i)
  {
    previous_i = current_i;
    vector<int>& current_neighs = VV[current_i];
    
    for(size_t i=0;i<current_neighs.size();++i)
    {
      double newv = squaredNorm(&(a_V.ptr()[current_neighs[i]*a_V_cols]),a_p.ptr());
      if (newv < current_v)
      {
        current_v = newv;
        current_i = current_neighs[i];
      }
    }
  }
  a_p.dealloc();
  return current_i;
}
int Phong::findClosestFace(const Vector8& p, int fid)
{
  igl::aligned_matrix a_p((RowVector8)(3*p.transpose()));
  
  RowVector8 target;
  igl::aligned_matrix a_target;
  
  target = (_V.row(_F(fid,0)) + _V.row(_F(fid,1)) + _V.row(_F(fid,2)) );
  a_target = igl::aligned_matrix(target);
  
  float current_v = squaredNorm(a_target.ptr(),a_p.ptr());
  int current_i = fid;
  int previous_i = -1;
  
  while (current_i != previous_i)
  {
    previous_i = current_i;
    const Eigen::Matrix<int, 1, 3> & current_neighs = TT.row(current_i);
    
    for(int i=0;i<3;++i)
    {
      if(current_neighs[i] == -1)
        continue;
      target = (_V.row(_F(current_neighs[i],0)) + _V.row(_F(current_neighs[i],1)) + _V.row(_F(current_neighs[i],2)) );
      a_target = igl::aligned_matrix(target);
      double newv = squaredNorm(a_target.ptr(),a_p.ptr());
      if (newv < current_v)
      {
        current_v = newv;
        current_i = current_neighs[i];
      }
    }
  }
  a_p.dealloc();
  a_target.dealloc();
  return current_i;
}

#ifdef ICC
float Phong::squaredNorm(const float* Erow, const float* p)
{
	__m256 avx_Erow, avx_p;
  
	// Erow-p
	avx_Erow = _mm256_load_ps(Erow);
	avx_p = _mm256_load_ps(p);
	avx_Erow = _mm256_sub_ps(avx_Erow,avx_p);
	
	// Erow-p * Erow-p
	avx_Erow = _mm256_dp_ps(avx_Erow,avx_Erow,0xF1);
	float f[8];
	_mm256_store_ps(f,avx_Erow);
	
	return f[0]+f[4];
}
#else
float Phong::squaredNorm(const float* Erow, const float* p)
{
  float acc = 0;
  float temp;
  // Erow - p * Erow - p
  for(int i=0; i<8;++i)
  {
    temp = Erow[i] - p[i];
    acc += temp * temp;
  }
  return acc;
}
#endif



void Phong::precomputeKRings()
{
  int K = 8;
  
    KRings.clear();
    KRings.resize(_F.rows());
    
    // Compute k-rings of vertices
    std::tr1::unordered_set<unsigned>myneighbors;
    std::pair<std::tr1::unordered_set<unsigned>::iterator, bool> inserted;
    std::tr1::unordered_set<unsigned>::iterator it;
    std::vector<unsigned> tocheck;
    // Initially every set contains only itself
    for (int i=0; i<_F.rows(); i++)
    {
      myneighbors.clear();
      KRings[i].clear();
      tocheck.clear();
      
      myneighbors.insert(i);
      KRings[i].push_back(i);
      tocheck.push_back(i);
      
      //dilate K times 
      for (int k=1; k<K; ++k)
      {
        vector<unsigned> prev = tocheck;
        tocheck.clear();
        for (unsigned int j = 0; j < prev.size(); ++j)
        {
          // add all the neighbours
          int fid = prev[j];
          for (int l=0;l<3;++l)
          {
            int toadd = TT(fid,l);
            if(toadd != -1)
            {
              inserted = myneighbors.insert(toadd);
              if(inserted.second)
              {
                tocheck.push_back(toadd);
                KRings[i].push_back(toadd);
              }
            }
          }
        }
      }
    }
}
