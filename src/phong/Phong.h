#ifndef PHONG_H
#define PHONG_H

#ifdef NO_LIBRARY_EXPORTS
#    define LIBRARY_API
#else
#    define LIBRARY_API __declspec(dllexport)
#endif

#include <assert.h>
#include "core/types.h"
#include <vector>
#include "aligned_matrix.h"
#include <chrono>

#ifdef ICC
#include <immintrin.h>
#endif

#include "core/PointOnSurface.h"
#include "phong/TrianglePhong.h"

#define EPSNEG -0.0001

class Anchor_Set_Phong
{
public:
  Anchor_Set_Phong(const std::vector<PointOnSurface > &point_set, const MatrixXX& E, const MatrixX3i& F);
  ~Anchor_Set_Phong();
  
  void dealloc(){a_anchors_E.dealloc();};
  
  MatrixXX anchors_E;
  igl::aligned_matrix a_anchors_E;
  std::vector<PointOnSurface > _point_set;
  
};

class Phong
{
public:
  typedef Eigen::Matrix<double, 2, 8> Basis;
  typedef Eigen::Matrix<ScalarType, 8, 1> Vector8;
  typedef Eigen::Matrix<IndexType, 1, Eigen::Dynamic> RowVectorXi;

  
  Phong();
  
  ~Phong();

  void init(const MatrixXX& V, const MatrixXXi& F, std::vector<Basis>* basis = 0, bool x_bUseEuclideanInit = false);
  
  // Project a point in the ambient space onto a triangle and return
  // the barycentric coordinates of the projection
  // Input:
  // fid is the id of the triangle
  // p is the point to project
  bool projectOnTriangle(const unsigned &fid, const Vector8 &p, RowVector3& w);
  
  // Project a point in ambient space to the mesh
  bool projectBruteForce(const Vector8 &p, int& fid, RowVector3& w, bool Phong);
  
  bool project(const Vector8 &p, int fid_start, int& fid, RowVector3& w);
  bool projectOnTriangle_fast(const unsigned &fid,
                              const float *p,
                              float *w);
  bool projectOnTriangleEuclidean_fast(const unsigned &fid,
                              const float *p,
                              float *w);
  
  // Unproject a point, providing barycentric coordinates for a point p that projects inside the provided point inside a triangle
  void unproject_fast(const Anchor_Set_Phong& anchor_set,
                      const PointOnSurface &point,
                      RowVectorX &W,
                      const double &ilya_hack_parameter = 10.);
  
  // Find the closest point on the surface starting from vid
  int findClosest(const Vector8& p, int vid);
  
  int findClosestFace(const Vector8& p, int vid);
  
  // find optimum by jumping across faces using the least negative weight
  bool findMinimumJumpNeg(const float *p, int fid_start, int& fid, float *w);
  
  bool projectBFS(const float *p, int fid_start, int& fid, float* w);
  
  std::vector<igl::aligned_matrix> _a_basis;
  Eigen::Matrix<Eigen::Matrix<ScalarType,2,9>, Eigen::Dynamic, 3> _baryop_cropped;
  // Barycentric operator
  igl::aligned_matrix _a_baryop_col_07;
  igl::aligned_matrix _a_baryop_col_8;
  igl::aligned_matrix _a_baryop_col_8_f;
  // Evaluate averages of points in E, with the weights in v
private:
  
  // -------------- Types

  // Type for the base of the surface tangent space at vertices

  // -------------- Functions
  
  // Compute a basis for the tangent space at every vertex
  std::vector<Basis> computeBasis(const MatrixXX& V, const MatrixXXi& F);

  std::vector<Basis> computeBasisE(const MatrixXX& V, const MatrixXXi& F);
  
  // Precomputation of barycentric projection operator
  void precomputeBarycentricOperator();
  
  // Precompute K-rings
  void precomputeKRings();
  
  // Compute squared norm of (Erow-p)
  float squaredNorm(const float* Erow, const float* p);
  
#ifdef ICC
  inline __m128 CrossProduct(__m128 a, __m128 b)
  {
    return _mm_sub_ps(
                      _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 1, 0, 2))), 
                      _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1)))
                      );
  };

  inline __m256 TwoCrossProducts(__m256 a, __m256 b)
  {
    return _mm256_sub_ps(
                         _mm256_mul_ps(_mm256_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1)), _mm256_shuffle_ps(b, b, _MM_SHUFFLE(3, 1, 0, 2))), 
                         _mm256_mul_ps(_mm256_shuffle_ps(a, a, _MM_SHUFFLE(3, 1, 0, 2)), _mm256_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1)))
                         );
  };

//  inline __m256 TwoCrossProducts(__m256 a, __m256 b)
//  {
//    return _mm256_sub_ps(
//                         _mm256_mul_ps(_mm256_shuffle_ps(a, a, _MM_SHUFFLE(3, 1, 0, 2)), _mm256_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1))), 
//                         _mm256_mul_ps(_mm256_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1)), _mm256_shuffle_ps(b, b, _MM_SHUFFLE(3, 1, 0, 2)))
//                         );
//  };
#else
  // evaluate cross product
  inline void cross(const float *a, const float *b, float *c)
  {
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
  };
  
  inline void crossd(const double *a, const double *b, double *c)
  {
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
  };
#endif
  // -------------- Fields

  // A base for every vertex
  std::vector<Basis> _basis;
  
  // Barycentric operator for vertices
  Eigen::Matrix<Eigen::Matrix<ScalarType,3,9>, Eigen::Dynamic, 3> _baryop;
  
  // Mesh representation
  MatrixXX  _V;
  MatrixXXi  _F;
  
  // AVX Aligned matrices
  igl::aligned_matrix a_V;
  unsigned a_V_cols;
  
  // VV relation
  std::vector<std::vector<int> > VV;
  
  // TT relation
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>  TT, TTi;

  // VF relation
  std::vector<std::vector<int> > VF;
  std::vector<std::vector<int> > VFi;
  
  // Sorted k-rings
  std::vector<std::vector<int> > KRings;
  
public:
  TrianglePhong TrianglesPhong;
  double initTime;

  int sizeV()
  {
	  return _V.rows();
  }

  int sizeF()
  {
	  return _F.rows();
  }
    
};

extern "C"
{
	LIBRARY_API Phong* createPhongObject(double* V, const int nV, const int dim, unsigned int* F, const int nF);
	LIBRARY_API bool project(Phong *phong, const double* p, int fid_start, float *w);
	LIBRARY_API bool deletePhongObject(Phong *phong);
	LIBRARY_API double getInitTime(Phong *phong);
}
#endif
