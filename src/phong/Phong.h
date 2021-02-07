#ifndef PHONG_H
#define PHONG_H

#ifdef NO_LIBRARY_EXPORTS
#    define LIBRARY_API
#else
	#ifdef _MSC_VER
		#define LIBRARY_API __declspec(dllexport)
	#else
		#define LIBRARY_API __attribute__((visibility("default")))
	#endif
#endif

#include <assert.h>
#include "core/types.h"
#include <vector>
#include "aligned_matrix.h"
#include <chrono>

#ifdef __INTEL_COMPILER
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
  
#ifdef __INTEL_COMPILER
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
	/*
	Initialize Phong projection for a triangle mesh (V, F)
	V: vertices [v11, ... v1n, v21, ... v2n, ..., v81, ... v8n]
	nV: number of vertices (n above)
	dim: dimension of Euclidean space in which V is embedded (must be 8)
	F: faces [F11, ..., F1m, F21, ..., F2m, F31, ..., F3m]
	nF: number of faces (m above)
	*/
	LIBRARY_API Phong* createPhongObject(double* V, const int nV, const int dim, unsigned int* F, const int nF);
	
	/*
	Project a point onto the triangle mesh. Returns the implicit coordinate (triangle index, barycentric coordinate) of the projection.
	phong: Pointer to a Phong object (typically returned by createPhongObject)
	p: point in dim-dimensional space [p1, ... p8]
	fid_start: index of face to start from.
	w: output variable: w[0]: face_index, (w[1], w[2], w[3]): barycentric 
	coordinate.
	*/
	LIBRARY_API bool project(Phong *phong, const double* p, int fid_start, float *w);
	
	/*
	Project a point onto the triangle mesh. Returns the implicit coordinate (triangle index, barycentric coordinate) of the projection.
	project() uses a short cut: starting from fid_start, find the first triangle
	in breadth-first order which returns a valid Phong projection for the input
	point. In contract, phongBruteForce() computes the Phong projection for each
	triangle, and returns the valid project that is nearest to the input point.
	projectBruteForce() can be more accurate than project() but slower.
	phong: Pointer to a Phong object (typically returned by createPhongObject)
	p: point in dim-dimensional space [p1, ... p8]
	w: output variable: w[0]: face_index, (w[1], w[2], w[3]): barycentric
	coordinate.
	*/
	LIBRARY_API bool projectBruteForce(Phong* phong, const double *p, float *w);
	
	/* 
	Returns the Euclidean (closest-point) projection on the mesh.
	phong: Pointer to a Phong object (typically returned by createPhongObject)
	p: point in dim-dimensional space [p1, ... p8]
	w: output variable: w[0]: face_index, (w[1], w[2], w[3]): barycentric coordinate.
	*/
	LIBRARY_API bool projectBruteForceEuclidean(Phong* phong, const double *p, float *w);
	
	/* 
	Deletes the Phong object.
	phong: Pointer to a Phong object (typically returned by createPhongObject)
	*/
	LIBRARY_API bool deletePhongObject(Phong *phong);
	
	/* 
	Returns the amount of time required to create the Phong object (in seconds).
	phong: Pointer to a Phong object (typically returned by createPhongObject)
	*/
	LIBRARY_API double getInitTime(Phong *phong);
}
#endif
