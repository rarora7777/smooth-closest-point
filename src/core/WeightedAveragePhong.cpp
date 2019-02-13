#include <iostream>

#include "WeightedAveragePhong.h"
#include "igl/mvc.h"
#include "igl/grad.h"
#include "igl/moveFV.h"

#ifdef ICC
#include <immintrin.h>  
#endif

#include "igl/Timer.h"


WeightedAveragePhong::WeightedAveragePhong(){};

WeightedAveragePhong::~WeightedAveragePhong(){};

WeightedAveragePhong::WeightedAveragePhong(const Mesh* mesh_,const MatrixXX embedding_matrix):
m_mesh_ptr(mesh_),
has_valid_mesh(false),
horizontal_parameter(10.)
{
  has_valid_mesh = m_mesh_ptr->loaded();
  embedding.init(embedding_matrix);
  if(m_mesh_ptr->loaded())
    if (m_mesh_ptr->num_vertices())
      phong.init(embedding.embedded_coordinates_const(),m_mesh_ptr->faces_const(), NULL, false);
};

// --------------------- Forward

bool WeightedAveragePhong::forward(const Anchor_Set_Phong &anchor_set,
                                   const RowVectorX &weights,
                                   PointOnSurface &point,
                                   int warm_start)
{
  assert (embedding.numPoints() == m_mesh_ptr->num_vertices());
  assert (anchor_set.anchors_E.rows() == size_t(weights.cols()));
  
  int closest_fid = -1;
  RowVector3 closest_w;
  
  Phong::Vector8 p;
  
  igl::aligned_matrix a_weights(weights);
  
  float a_p[8];
  PointOnSurface::evaluate_average(a_weights.ptr(), weights.size(), anchor_set.a_anchors_E.ptr(), a_p);

  for (size_t i = 0; i< 8; ++i)
    p[i] = a_p[i];
  a_weights.dealloc();
  
  //Use closest anchor as starting point
  //WARM START IS A FACE
  int fid;
  if (warm_start >=0 )
    fid = warm_start; // WARM START
  else
  {
    int closest_anchor;
    weights.maxCoeff(&closest_anchor); // CLOSEST ANCHOR
    fid = anchor_set._point_set[closest_anchor].fid;
  }
  
  bool b;
  
    b = phong.project(p, fid, closest_fid, closest_w);

  if (!b)
    return false;
  
  point = PointOnSurface(closest_fid, closest_w);

  
  return true;
};


// ---------------------- Backward

bool WeightedAveragePhong::backward(const Anchor_Set_Phong& anchor_set,
                                    const PointOnSurface &point,
                                    RowVectorX &weights)
{
  
  phong.unproject_fast(anchor_set,
                       point,
                       weights);
  return true;
}


// ------------------- Utilities

const Embedding& WeightedAveragePhong::embedding_const() const
{
  return embedding;
}