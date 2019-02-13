#ifndef Preview3D_lscoord_h
#define Preview3D_lscoord_h

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include "../../igl_lib/min_quad_with_fixed.h"

int lscoord( const Eigen::RowVectorXd &q, const Eigen::MatrixXd &P, Eigen::MatrixXd &W, const bool positive = false, const bool weighted = true)
{
  Eigen::MatrixXd G (P.rows(), P.cols());
  for (int i = 0; i<P.rows(); ++i)
    G.row(i) = P.row(i) - q;
  
  // Minimize |w|^2
  int n = P.rows();
  int m = P.cols();
  
  Eigen::SparseMatrix<double, Eigen::RowMajor>H = Eigen::SparseMatrix<double, Eigen::RowMajor> (n,n);
  H.reserve (n);
  for (int i =0 ; i<n; ++i)
  {
    H.startVec(i);
    if (weighted)
      H.insertBack(i,i) = (G.row(i)).dot(G.row(i));
    else
      H.insertBack(i,i) = 1.;
  }
  H.finalize();
  
  Eigen::MatrixXd f = Eigen::MatrixXd::Zero(n,1);
  
  
  Eigen::SparseMatrix<double, Eigen::RowMajor>Aeq = Eigen::SparseMatrix<double, Eigen::RowMajor> (m+1,n);
  Aeq.reserve (n*(m+1));
  Eigen::MatrixXd beq = Eigen::MatrixXd::Zero(m+1,1);
  
  // Must sum to 1...
  Aeq.startVec(0);
  for (int i = 0; i < n ; ++i)
    Aeq.insertBack(0,i) = 1.; 
  beq(0,0) = 1.;
  // sum_i w_i G = 0
  for (int j = 0; j < m ; ++j)
  {
    Aeq.startVec(j+1);
    for (int i = 0; i < n ; ++i)
      Aeq.insertBack(j+1,i) = G(i,j); 
    beq(j+1,0) = 0.;
  }
  Aeq.finalize();
  
  if (positive)
  {
    // Enforce positivity
    Eigen::MatrixXd A = -1*Eigen::MatrixXd::Identity(n,n);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
    
    // QP Solve
    //W = quadprog(H,f,A,b,Aeq,beq);
  }
  else
  {
    igl::min_quad_with_fixed_data <double> data;
    igl::min_quad_with_fixed_precompute <double> ( H, Eigen::MatrixXi(0,1), Aeq,  false, data );
    igl::min_quad_with_fixed_solve <double> ( data, f, Eigen::MatrixXd(0,1), beq, W);
    
  }
  return int((Aeq*W -beq).array().abs().matrix().minCoeff() >1e-10);
}


#endif
