#ifndef PointOnSurface_h
#define PointOnSurface_h

#include "core/types.h"
#include <vector>

using namespace std;

#include "phong/Phong.h"
class PointOnSurface
{
  friend struct PointOnSurfaceFinder;
public:
  PointOnSurface():
  fid(-1),
  has_e(false)
  {};
  
  PointOnSurface(const IndexType fid_, const Vector3 &bc_):
  bc(bc_/bc_.sum()),
  fid(fid_),
  has_e(false)
  {
    bc_ptr[0] = bc[0];
    bc_ptr[1] = bc[1];
    bc_ptr[2] = bc[2];
  };
  
  PointOnSurface(IndexType vid_, const MatrixX3i& F, const vector<vector< int > >& VF, const vector<vector< int > >& VFi) :
  has_e(false)  
  {
    fid = VF[vid_][0];
    bc = Vector3::Zero();
    assert(F(VF[vid_][0],VFi[vid_][0]) == vid_);
    bc(VFi[vid_][0]) = 1.0;
    
    bc_ptr[0] = bc[0];
    bc_ptr[1] = bc[1];
    bc_ptr[2] = bc[2];
  };
  
  ~PointOnSurface()
  {
  }
  // Variables
  
  VectorX get_position(const MatrixXX& V, const MatrixX3i& F) const
  {
    if (fid == -1)
      return VectorX::Zero(V.cols());
    else
      return 
      (
       V.row(F(fid,0)) * bc(0) +
       V.row(F(fid,1)) * bc(1) +
       V.row(F(fid,2)) * bc(2)
       );
  }
  
  void precompute(const MatrixXX& E, const MatrixX3i& F)
  {
    igl::aligned_matrix a_v((MatrixXX(3,8) << E.row(F(fid,0)),
                                              E.row(F(fid,1)),
                                              E.row(F(fid,2)) ).finished());
    evaluate_average(bc_ptr, 3, a_v.ptr(), e);
    a_v.dealloc();
    has_e = true;
  }
  
  int fid;
  Vector3 bc;
  float bc_ptr[3];
  
  bool has_e;

  float e[8];
  
#ifdef ICC
  static void evaluate_average(const float* v,
                        size_t v_cols,
                        const float* E,
                        float* accret)
  {
    __m256 E_row, v_row, temp;
    __m256 acc = _mm256_setzero_ps();
    
    // v*E
    for(int i=0; i<v_cols;++i)
    {
      E_row  = _mm256_load_ps(&(E[i * 8]));
      v_row  = _mm256_set1_ps(v[i]);
      temp   = _mm256_mul_ps(E_row,v_row);
      acc    = _mm256_add_ps(temp,acc);
    }
    
    _mm256_store_ps(accret,acc);
  }
#else
  static void evaluate_average(const float* v,
                        size_t v_cols,
                        const float* E,
                        float* acc)
  {
    for(int i=0; i<8;++i)
      acc[i] = 0;
    
    // v*E
    float temp;
    for(int a=0;a<8;++a)
    {
      temp = 0;
      for(size_t i=0; i<v_cols;++i)
      {
        temp   = E[i * 8 + a] * v[i];
        acc[a] += temp;
      }
    }
    
  }
#endif
  
  int findInVector (const std::vector<PointOnSurface> &points)
  {
    for (int i = 0; i<points.size(); ++i)
      if (points[i].fid == fid && (points[i].bc - bc).cwiseAbs().maxCoeff() <1e-4)
        return i;
    return -1;
  }
  
};


#endif
