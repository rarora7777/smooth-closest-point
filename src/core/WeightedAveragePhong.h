#ifndef Preview3D_WeightedAveragePhong_h
#define Preview3D_WeightedAveragePhong_h
#include "Mesh.h"
#include "Embedding.h"
#include "phong/Phong.h"


class WeightedAveragePhong
{
  friend class Relaxation_App;
  
public:
  // -------------- Constructors
  WeightedAveragePhong();
  
  ~WeightedAveragePhong();
  
  WeightedAveragePhong(const Mesh* mesh_,const MatrixXX embedding);
  
  // -------------- Forward
  
  bool forward(const Anchor_Set_Phong &anchor_set,
               const RowVectorX &weights_,
               PointOnSurface &point,
               int warm_start = -1);
  
  // --------------- Backward
  Anchor_Set_Phong generateAnchorSet(const std::vector<PointOnSurface > &anchor_set);
  
  bool backward(const Anchor_Set_Phong& anchor_set,
                const PointOnSurface &point,
                RowVectorX &weights);
    
  
  // -------------- Utilities
  void reset_use_brute_force_search(bool value);
  void reset_use_mvc(bool value);
  void reset_use_phong(bool value);
  void reset_use_gapfiller(bool value);
  const Embedding& embedding_const() const;
  
  const MatrixXX& Ec() const
  { return embedding.embedded_coordinates_const(); }
  const MatrixX3& Vc() const
  { return m_mesh_ptr->vertices_const(); }
  const MatrixX3i& Fc() const
  { return m_mesh_ptr->faces_const(); }
  const vector<vector<int> >& VFc() const
  { return m_mesh_ptr->VF; }
  const vector<vector<int> >& VFic() const
  { return m_mesh_ptr->VFi; }

  
  
  //private:
  const Mesh *m_mesh_ptr;
  Embedding embedding;
  
  bool has_valid_mesh;
  
  double horizontal_parameter;
  
  Phong phong;
    
};


#endif
