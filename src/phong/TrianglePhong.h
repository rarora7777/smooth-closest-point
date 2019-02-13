/*--
 TrianglePhong.h
 
 All matrices are stored in Row-Major order
 */

#ifndef TRIANGLEPHONG_H_INCLUDED
#define TRIANGLEPHONG_H_INCLUDED

#include <Eigen/Core>
#include "aligned_matrix.h"

#define N_ITER 10
//#define TIMERS_ENABLED

class TrianglePhong
{
public:
  typedef Eigen::Matrix<float, 8, 1> Vector;
  typedef Eigen::Vector3f Vector3;
  typedef Eigen::Matrix2f Matrix2;
  typedef Eigen::Matrix3f Matrix3;
  typedef Eigen::RowVector3f RowVector3;
  typedef Eigen::Matrix<float, 2, 8> Basis;
  
  TrianglePhong() {}
  // Preprocess a mesh for Phong Projection
  // The mesh must be embedded in a 8d space, if you want to use this class for a 3d mesh just pad with zeros
  // Corners: #F x 24 matrix with the coordinates of the vertices of every triangle in every row
  //          in the format x1, y1, z1, ...,  x2, y2, z2, ..., x3, y3, z3, ...
  // Tangents: #F x 48 matrix with the tangent planes, every row is associated to a face, and there are three
  //           tangent planes for every face, each of them composed of two vectors of 8 elements
  //           Every row has the format: tp1row1 tp1row2 tp2row1 tp2row2 tp3row1 tp3row2
  TrianglePhong(Eigen::MatrixXd Corners, Eigen::MatrixXd Tangents);

  // Projects a point onto the triangle fid, point must be a float[8] and w a float[3]
  // If the projection is invalid, it returns -1 -1 3
  void project(unsigned fid, const float* point, float* w);
  
  // Blend the tangent planes of triangle fid with weights bary
  // bary is a float[3], blended a float[16]
  void blendf(unsigned fid, const float* bary, float* blended);
  
  // Blend the vertices position at triangle fid with weights bary
  // bary is a float[3], blended a float[8]
  void blendPosf(unsigned fid, float* bary, float* blended);

  
private:
  
  void blend_aux(unsigned fid, float coeff[6], float* blended);
  void computeW(const float* bary, float coeff[6]);
  void precompute(Vector corners[3], Basis tangents[3], Basis* ert, Basis* et1);
  
  // #V x 24: mesh vertices per triangle
  igl::aligned_matrix _corners;
  float* getCorner(unsigned fid,unsigned i)
  {
    assert(i>=0 && i<=3);
    return _corners.ptr() + fid*24 + i*8;
  }

  // #V x 48: 3 tangent planes 2x8 per triangle
  igl::aligned_matrix _tangents;
  float* getTangent(unsigned fid,unsigned i)
  {
    assert(i>=0 && i<=3);
    return _tangents.ptr() + fid*48 + i*16;
  }

  //Matrix2 _r[3]; //rotations such that r[i] * t[i] is similar to t[i + 1]
  //Matrix2 _e[3]; //rotations such that e[i] * r[i] * t[i] and e[i] * t[i + 1] are similar for all i.
  
  // _ert = (_e[i] * _r[i])).eval() * _tangents[i]
  // #V x 48: 3 tangent planes 2x8 per triangle
  igl::aligned_matrix _ert;
  float* getErt(unsigned fid,unsigned i)
  {
    assert(i>=0 && i<=3);
    return _ert.ptr() + fid*48 + i*16;
  }

  // _et1 = _e[i] * _tangents[(i+1)%3]
  // #V x 48: 3 tangent planes 2x8 per triangle
  igl::aligned_matrix _et1;
  float* getEt1(unsigned fid,unsigned i)
  {
    assert(i>=0 && i<=3);
    return _et1.ptr() + fid*48 + i*16;
  }

  // Precomputed stuff for first Newton Iteration
  igl::aligned_matrix _blendB;
  igl::aligned_matrix _blendp;
  igl::aligned_matrix _minusbv;
  igl::aligned_matrix _blended;
  
};


#endif //TRIANGLEPHONG_H_INCLUDED
