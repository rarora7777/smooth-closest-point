#ifndef Preview3D_Mesh_h
#define Preview3D_Mesh_h

#include <vector>
#include "types.h"

#include "igl/adjacency_list.h"
#include "igl/vf.h"

#include <iostream>

using namespace std;

class Mesh
{
  friend class Relaxation_App;
  
public:
  Mesh(){is_loaded = false;};
  
  Mesh(const PointMatrixType* vertices_,
       const FaceMatrixType* faces_):
  face_size(0),
  is_loaded(false)
  {
    if(vertices_->cols() != 3)
    {
      std::cerr<< "Mesh() : Error reading vertices: Each vertex should be a 3-vector."<<std::endl;
      exit(-1);
    }
    vertices = *vertices_;
    
    if(faces_->rows())
    {
      face_size = faces_->cols();
      if(face_size !=3 && face_size !=4)
      {
        std::cerr<< "Mesh() : Error reading mesh: Each face should be a triangle or a quadrangle."<<std::endl;
        exit(-1);
      }
    }    
    
    if(face_size == 3)
      faces = *faces_;
    else
    {
      faces = MatrixXXi(2*faces_->rows(),3);
      faces.block(0,0,faces_->rows(), 3) = *faces_;
      faces.block(faces_->rows()+1,0,faces_->rows(), 1) = faces_->col(0);
      faces.block(faces_->rows()+1,1,faces_->rows(), 1) = faces_->col(2);
      faces.block(faces_->rows()+1,2,faces_->rows(), 1) = faces_->col(3);
      quad_faces = *faces_;
    }    
    
    // Init

    if (vertices.rows() > 0)
    {
      compute_bounding_box_and_centroid();
      compute_vertex_neighbors();
      compute_face_normals();
      compute_vertex_normals();
      
      get_min_max_edge();
      get_avg_edge();
    }
    
    // End Init
    
    is_loaded = vertices.rows()>0;
  }
  ~Mesh(){};
  
  template<typename Mat>
  void return_vertex_positions(Mat *vertices_v)
  {
    if(!vertices_v)
      return;
    vertices_v->resize(vertices.rows(),3);
    for(int i =0 ; i<vertices.rows(); ++i)
    {
      vertices_v->row(i) <<vertices(i,0), vertices(i,1), vertices(i,2);
    }
  }  
  
  bool loaded() const 
  {
    return is_loaded;
  }
  
  double min_edge() const
  {
    return min_edge_length;
  }
  
  double max_edge() const
  {
    return max_edge_length;
  }
  
  double avg_edge() const
  {
    return avg_edge_length;
  }
  
  double min_radius() const 
  {
    return 0.5*(bounding_box.row(1) - bounding_box.row(0)).cwiseAbs().minCoeff();
  }
  
  double max_radius() const 
  {
    return 0.5*(bounding_box.row(1) - bounding_box.row(0)).cwiseAbs().maxCoeff();
  }
  
  double diameter() const
  {
    return (bounding_box.row(1) - bounding_box.row(0)).norm();
  }
  
  int num_vertices() const {return vertices.rows();};
  int num_faces() const {return faces.rows();};
    
  const std::vector<int>& vertex_neighbors(const int id) const
  {
    assert (id>=0 && id<vertices.rows());
    return VV[id];
  }
  
  const std::vector<int>& face_neighbors(const int id) const
  {
    assert (id>=0 && id<vertices.rows());
    return VF[id];
  }
  
  const MatrixX3& vertices_const() const
  {
    return vertices;
  }
  
  const MatrixX3i& faces_const() const
  {
    return faces;
  }
  
  //private:
  
  MatrixX3 vertices;
  MatrixX3 vertex_normals;
  
  //only triangles allowed
  MatrixX3i faces;
  MatrixX4i quad_faces;
  MatrixX3i face_texture_indices;
  MatrixX3 face_normals;
  
  int face_size;
  
  std::vector<std::vector<int> > VF;
  std::vector<std::vector<int> > VFi;
  
  std::vector<std::vector<int> > VV;
  
  RowVector3 centroid;
  Eigen::Matrix<ScalarType, 2, 3> bounding_box;
  double min_edge_length, max_edge_length, avg_edge_length;
  
  bool is_loaded;
  
  void compute_face_normals()
  {
//    printf("Calculating Face Normals...\n");
    face_normals = MatrixXX(faces.rows(), 3);
    // loop over faces
    for(int i = 0; i <faces.rows(); ++i)
    {
      // Only valid for triangle meshes !!!
      RowVector3 v1 = vertices.row(faces(i,1)) - vertices.row(faces(i,0));
      RowVector3 v2 = vertices.row(faces(i,2)) - vertices.row(faces(i,0));
      face_normals.row(i) = v1.cross(v2).normalized();
    }
  }
  
  void compute_vertex_normals()
  {
//    printf("Calculating Vertex Normals...\n");
    //    std::vector<Vector > face_normals;
    if (face_normals.rows() == 0)
      compute_face_normals();
    
    vertex_normals = MatrixXX::Zero(vertices.rows(), 3);
    if(false)
    {
      
      // loop over faces
      for(int i = 0; i <faces.rows(); ++i)
        // loop over vertices in this face
        for(int j = 0; j <faces.cols(); ++j)
          vertex_normals.row(faces(i,j)) += face_normals.row(i);
    }
    
    else
    {
      vector< vector<int> > VF;
      vector< vector<int> > VFi;  
      
      igl::vf(vertices,faces,VF,VFi);
      for(int i = 0; i <vertices.rows(); ++i)
      {
        const RowVectorX &p = vertices.row(i);
        for(unsigned j = 0; j <VF[i].size(); ++j)
        {
          const RowVector3i &face = faces.row(VF[i][j]);
          int next = (VFi[i][j]+1)%3; const RowVectorX &next_vertex = vertices.row(face[next]);
          int prev = (VFi[i][j]+2)%3; const RowVectorX &prev_vertex = vertices.row(face[prev]);
          double angle = std::acos((next_vertex - p).dot(prev_vertex - p));
          double w = angle / (next_vertex - p).norm() / (prev_vertex - p).norm();
          vertex_normals.row(i) += w*face_normals.row(VF[i][j]);
        }
      }
    }

    for(int i = 0; i <vertices.rows(); ++i)
      vertex_normals.row(i).normalize();
  }
  
  //vertex_neighbors, unsorted
  void compute_vertex_neighbors()
  {
    igl::vf(vertices,faces,VF,VFi);
    igl::adjacency_list(faces, VV,true);
  }
  
  void compute_bounding_box_and_centroid()
  {
    
    centroid = RowVector3::Zero();
    bounding_box.row(0) = RowVector3::Constant(std::numeric_limits<ScalarType>::max());
    bounding_box.row(1) = RowVector3::Constant(-std::numeric_limits<ScalarType>::max());
    for(int vi = 0; vi < vertices.rows(); ++vi)
    {
      centroid += vertices.row(vi);
      bounding_box.row(0) = bounding_box.row(0).array().min(vertices.row(vi).array()).matrix(); 
      bounding_box.row(1) = bounding_box.row(1).array().max(vertices.row(vi).array()).matrix();
    }
    centroid /= vertices.rows();
  }
  
  void get_min_max_edge()
  {
    min_edge_length = 1e10;
    max_edge_length = 0;
    for (int i=0; i<vertices.rows(); ++i)
    {
      for (unsigned j = 0; j<VV[i].size(); ++j)
      {
        min_edge_length = std::min<ScalarType>( (vertices.row(i) - vertices.row(VV[i][j])).squaredNorm(), min_edge_length );
        max_edge_length = std::max<ScalarType>( (vertices.row(i) - vertices.row(VV[i][j])).squaredNorm(), max_edge_length );
      }
    }
    min_edge_length = std::sqrt(min_edge_length);
    max_edge_length = std::sqrt(max_edge_length);
  }
  
  void get_avg_edge()
  {
    avg_edge_length = 0;
    double num = 0;
    for (unsigned i=0; i<vertices.rows(); ++i)
    {
      for (unsigned j = 0; j<VV[i].size(); ++j, ++num)
      {
        double temp = (vertices.row(i) - vertices.row(VV[i][j])).norm();
        avg_edge_length += temp;
      }
    }
    avg_edge_length /= num;
  }
  
};


#endif
