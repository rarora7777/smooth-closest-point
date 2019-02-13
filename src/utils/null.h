#ifndef Preview3D_null_h
#define Preview3D_null_h

#include <Eigen/Core>

template <typename T>
inline void nullspace(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
                      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &Z)

{
  Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > svd(A, Eigen::ComputeFullV );
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> V = svd.matrixV();
  const Eigen::Matrix<T, Eigen::Dynamic, 1> s = svd.singularValues();
  
  T tol = max(A.rows(),A.cols()) * s.maxCoeff() * std::numeric_limits<T>::epsilon();
  int r;
  for (r = 0; r < s.rows(); ++r)
  {
    if (s[r] < tol)
      break;
  }

  //keep r last columns of V
  Z = (V.block(0,r,V.rows(),V.cols()-r));
}

#endif
