#ifndef Preview3D_Embedding_h
#define Preview3D_Embedding_h

#include "types.h"
#include "Mesh.h"
#include "aligned_matrix.h"


class Embedding
{
public:
  Embedding(){};
  ~Embedding(){};

  void init(const MatrixXX& mat)
  {
    assert(mat.cols() == 8);
    embedded_coordinates = mat;

    a_embedded_coordinates = igl::aligned_matrix(embedded_coordinates);
  };
  
  int dimension() const 
  {
    return embedded_coordinates.cols();
  }
  
  int numPoints() const 
  {
    return embedded_coordinates.rows();
  }
  
  const MatrixXX& embedded_coordinates_const() const
  {
    return embedded_coordinates;
  }
  
  //this used to be E in the Matlab code
  MatrixXX embedded_coordinates;
  
  igl::aligned_matrix a_embedded_coordinates;
  
};


#endif
