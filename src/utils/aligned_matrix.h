#ifndef IGL_ALIGNED_MATRIX_H
#define IGL_ALIGNED_MATRIX_H

#define posix_memalign(p, a, s) (((*(p)) = _aligned_malloc((s), (a))), *(p) ?0 :errno)

#include <vector>

namespace igl 
{
  class aligned_matrix
  {
  public:

    aligned_matrix()
    {
      _M = 0;
      _r = 0;
      _c = 0;
    }

    aligned_matrix(unsigned m, unsigned n)
    {
      _r = m;
      _c = n;
      
      _M = 0;
      
      if (posix_memalign((void**)&_M,32,4*n*m) != 0)
        assert(0);
    }

    
    aligned_matrix(const Eigen::MatrixXd& M)
    {
      _r = M.rows();
      _c = M.cols();
      
      _M = 0;
      if (posix_memalign((void**)&_M,32,4*M.cols()*M.rows()) != 0)
        assert(0);
      for(unsigned i=0;i<M.rows(); ++i)
        for(unsigned j=0;j<M.cols();++j)
          _M[i*M.cols() + j] = M(i,j);

    }

    aligned_matrix(const Eigen::Matrix<double, 2, 8>& M)
    {
      _r = M.rows();
      _c = M.cols();
      
      _M = 0;
      if (posix_memalign((void**)&_M,32,4*M.cols()*M.rows()) != 0)
        assert(0);
      for(unsigned i=0;i<M.rows(); ++i)
        for(unsigned j=0;j<M.cols();++j)
          _M[i*M.cols() + j] = M(i,j);
      
    }
    
    aligned_matrix(const Eigen::Matrix<double, 1, Eigen::Dynamic>& M)
    {
      _r = 1;
      _c = M.size();
      
      _M = 0;
      if (posix_memalign((void**)&_M,32,4*_r*_c) != 0)
        assert(0);
      for(unsigned j=0;j<_c;++j)
        _M[j] = M(j);
      
    }

    aligned_matrix(const Eigen::Matrix<double, Eigen::Dynamic, 1>& M)
    {
      _r = M.size();
      _c = 1;
      
      _M = 0;
      if (posix_memalign((void**)&_M,32,4*_r*_c) != 0)
        assert(0);
      for(unsigned j=0;j<_r;++j)
        _M[j] = M(j);
      
    }
    

    aligned_matrix(const Eigen::Matrix<double, 1, 8>& M)
    {
      _r = 1;
      _c = 8;
      
      _M = 0;
      if (posix_memalign((void**)&_M,32,4*_r*_c) != 0)
        assert(0);
      for(unsigned j=0;j<_c;++j)
        _M[j] = M(j);
      
    }


    
    ~aligned_matrix()
    {
    }
    
    void dealloc()
    {
      if (_M)
        free(_M);
    }
    
    inline size_t& rows()
    {
      return _r;
    }

    inline float* row(size_t i)
    {
      return &_M[i*_c];
    }

    inline size_t& cols()
    {
      return _c;
    }

    inline float * const ptr() const 
    {
      return _M;
    }
    
    inline float& at(size_t i, size_t j)
    {
      assert(i*_c + j >= 0);
      assert(i*_c + j < (_c*_r));
      return _M[i*_c + j];
    }
    
    size_t _r,_c;
    float* _M;
  };

  class aligned_row_list
  {
  public:
    
    aligned_row_list()
    {
      _r = 0;
    }
    
    aligned_row_list(Eigen::MatrixXd& M)
    {
      _r = M.rows();
      
      for(unsigned i=0; i<_r;++i)
      {
        // Allocate a single aligned row
        float* row;
        if (posix_memalign((void**)&row,32,4*M.cols()) != 0)
          assert(0);
        for(int j=0;j<M.cols();++j)
          row[j] = M(i,j);
        _c.push_back(M.cols());
        _M.push_back(row);
      }
      
    }
    
    ~aligned_row_list()
    {
//      if (_r != 0)
//        for(int i=0; i<_r; ++i)
//          free(_M[i]);
    }

    void dealloc()
    {
      for(unsigned i=0; i<_r; ++i)
        free(_M[i]);
    }
    
    size_t& rows()
    {
      return _r;
    }
    
    size_t& cols(size_t i)
    {
      return _c[i];
    }
    
    float* ptr(size_t i)
    {
      return _M[i];
    }
    
    float& at(size_t i, size_t j)
    {
      assert(j >= 0);
      assert(j < cols(i));
      return _M[i][j];
    }
    
    std::vector<float*> _M;
    std::vector<size_t> _c;
    size_t _r;
  };

}

#endif
