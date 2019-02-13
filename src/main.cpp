// Internal
#include "core/WeightedAveragePhong.h"
#include "core/Mesh.h"
#include "core/Embedding.h"

// System
#include <string>
#include <iostream>
#include <fstream>

// libIGL
#include <igl/readOBJ.h>

void help()
{
	cout << "Usage   : " << endl;
	cout << "Forward : WA f mesh.obj anchors weights out.txt" << endl;
	cout << "Inverse : WA i mesh.obj anchors points out.txt" << endl;
}

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

Eigen::MatrixXd readMatrix(string filename)
{
  // Read everything
  std::string line;
  std::ifstream infile(filename.c_str());
  
  vector<vector<double> > data;
  
  while (std::getline(infile, line))
  {
    std::istringstream iss(line);
    vector<double> data_line;
    double c;
    
    while (iss >> c)
      data_line.push_back(c);

    data.push_back(data_line);
  }
  infile.close();
  
  // Check if each row has the same number of entries
  int m = data.size();
  int n = data[0].size();
  
  for (unsigned i=0;i<data.size();++i)
    assert(data[i].size() == n);
  
  // Allocate an Eigen matrix
  Eigen::MatrixXd mat(m,n);
  for (unsigned i=0;i<m;i++)
    for (unsigned j=0;j<n;j++)
      mat(i,j) = data[i][j];
  
  return mat;
}

void writeMatrix(Eigen::MatrixXd mat, string filename)
{
  // Read everything
  std::ofstream outfile(filename.c_str());
  
  int m = mat.rows();
  int n = mat.cols();
  
  for (unsigned i=0;i<m;i++)
  {
    for (unsigned j=0;j<n;j++)
      outfile << mat(i,j) << " ";
    outfile << std::endl;
  }
  
  outfile.close();
}


int main(int argc, char **argv)
{
  if (argc != 6)
  {
    help();
    exit(0);
  }
	
  string mesh_path(argv[2]);
  string anchors_path(argv[3]);
  string weights_path(argv[4]);
  string out_path(argv[5]);
  Mesh M;
  
  // Load Mesh
  if (fexists(mesh_path.c_str()))
  {
    PointMatrixType V;
    FaceMatrixType F;
    igl::readOBJ(mesh_path, V, F);
    M = Mesh(&V,&F);
  }
  else
  {
    cout << "Mesh file missing or corrupted: " << mesh_path << endl;
    exit(1);
  }
  
  
  // Load Embedding
  Eigen::MatrixXd embedding;
  string embedding_path = mesh_path + string(".emb");
  if (fexists(embedding_path.c_str()))
  {
    embedding = readMatrix(embedding_path.c_str());
    
  }
  else
  {
    cout << "Embedding missing or corrupted: " << embedding_path << endl;
    exit(1);
  }
  
  // Load Anchors
  Eigen::MatrixXd anchors;
  if (fexists(anchors_path.c_str()))
  {
    anchors = readMatrix(anchors_path.c_str());
  }
  else
  {
    cout << "Anchors matrix missing or corrupted: " << anchors_path << endl;
    exit(1);
  }
  
  // Load Weights/Points
  Eigen::MatrixXd weights;
  if (fexists(weights_path.c_str()))
  {
    weights = readMatrix(weights_path.c_str());
  }
  else
  {
    cout << "Weights/Points matrix missing or corrupted: " << weights_path << endl;
    exit(1);
  }
  
  // Initialize WA
  WeightedAveragePhong WA(&M,embedding);

  // Initialize the anchors
  if (anchors.cols() != 3)
  {
    cout << "Each anchors should be specified as: fid bc1 bc2" <<  endl;
    exit(1);
  }

  vector<PointOnSurface> Anchors;
  for (unsigned i=0;i<anchors.rows();++i)
  {
    Vector3 v;
    v[0] = anchors(i,1);
    v[1] = anchors(i,2);
    v[2] = 1-(v[0] + v[1]);
    
    PointOnSurface p(anchors(i,0),v);
    Anchors.push_back(p);
  }
  
  Anchor_Set_Phong AS(Anchors,WA.Ec(),WA.Fc());
  
  if (argv[1][0] == 'f')
  {
    // Forward
    if (anchors.rows() != weights.cols())
    {
      cout << "Inconsistent input. Please provide as many weights as the number of anchors" << endl;
      exit(1);
    }

    vector<PointOnSurface> result;
    
    // Solve Forward
    for (unsigned i=0; i<weights.rows(); ++i)
    {
      PointOnSurface p;
      RowVectorX w = weights.row(i);
      WA.forward(AS,w,p);
      result.push_back(p);
    }

    // Export to a text file
    MatrixXX temp(result.size(),3);
    for (unsigned i=0; i<weights.rows(); ++i)
      temp.row(i) << result[i].fid, result[i].bc[0], result[i].bc[1];
    writeMatrix(temp, out_path);
  }
  else
  {
    // Inverse
    if (argv[1][0] != 'i')
    {
      help();
      exit(1);
    }
    
    if (weights.cols() != 3)
    {
      cout << "Inconsistent input. Please provide the points as: fid bc1 bc2" << endl;
      exit(1);
    }
    
    MatrixXX result(weights.rows(),anchors.rows());

    // Solve Inverse
    for (unsigned i=0; i<weights.rows(); ++i)
    {
      Vector3 v;
      v[0] = weights(i,1);
      v[1] = weights(i,2);
      v[2] = 1-(v[0] + v[1]);

      PointOnSurface p(weights(i,0),v);
      RowVectorX w;
      WA.backward(AS,p,w);

      result.row(i) = w;
    }
    
    // Export to a text file
    writeMatrix(result, out_path);
  }
}

