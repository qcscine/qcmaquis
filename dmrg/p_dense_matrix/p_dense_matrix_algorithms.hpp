#ifndef __ALPS_DENSE_MATRIX_ALGORITHMS_HPP__
#define __ALPS_DENSE_MATRIX_ALGORITHMS_HPP__

#include "p_dense_matrix/concept/matrix_concept_check.hpp"

#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include "p_dense_matrix/p_dense_matrix.h"
#include "utils/timings.h"

namespace blas
{
    template <typename Matrix>
    Matrix transpose(Matrix const& m) 
    {
        static Timer timer("transpose");
        timer.begin();
        BOOST_CONCEPT_ASSERT((blas::Matrix<Matrix>)); 
        // TODO: perhaps this could return a proxy object
        Matrix tmp(num_cols(m), num_rows(m));
        for(typename Matrix::size_type i=0; i < num_rows(m); ++i){
            for(typename Matrix::size_type j=0; j < num_cols(m); ++j){
                tmp(j,i) = m(i,j);
            }
        }
        timer.end();
        return tmp;
    }
        
    template<typename T>
    p_dense_matrix<T> conjugate(p_dense_matrix<T> M)
    {
        assert(false);
        M.inplace_conjugate();
        return M;
    }

    template <typename Matrix>
    const typename Matrix::value_type trace(Matrix const& m)
    {
        BOOST_CONCEPT_ASSERT((blas::Matrix<Matrix>)); 
        assert(num_rows(m) == num_cols(m));
        typename Matrix::value_type tr(m(0,0));
        for(typename Matrix::size_type i = 1; i<num_rows(m); ++i)
            tr += m(i,i);
        return tr;
    }
        
    template<class Matrix>
    Matrix identity_matrix(typename Matrix::size_type size)
    {
        assert(false);
        Matrix ret(size, size); // don't forget null-init here
        for (typename Matrix::size_type k = 0; k < size; ++k)
            ret(k,k) = 1;
        return ret;
    }
    
    template<typename T>
    void pblas_gemm(const p_dense_matrix<T>& A, const p_dense_matrix<T>& B, p_dense_matrix<T>& C)
    {
        C.resize(A.num_rows(), B.num_cols());
        C.set_init(ambient::null_i<T>);
	ambient::push(ambient::gemm_l_scalapack, ambient::gemm_c_scalapack, A, B, C);
    }

    template<typename T>
    void gemm(const p_dense_matrix<T>& A, const p_dense_matrix<T>& B, p_dense_matrix<T> &C)
    {
        C.resize(A.num_rows(), B.num_cols());
        C.set_init(ambient::null_i<T>);
	ambient::push(ambient::gemm_l, ambient::gemm_c, A, B, C);
        ambient::playout();
    }

    template<typename T>
    void gemm(const p_dense_matrix<T>& A, const p_diagonal_matrix<T>& B, p_dense_matrix<T>& C)
    {
        assert(num_cols(A) == num_rows(B));
        C.resize(A.num_rows(), B.num_cols());
	ambient::push(ambient::gemm_diagonal_rhs_l, ambient::gemm_diagonal_rhs_c, A, B.get_data(), C);
        ambient::playout();
    }
    
    template<typename T>
    void gemm(const p_diagonal_matrix<T>& A, const p_dense_matrix<T>& B, p_dense_matrix<T>& C)
    {
        assert(num_cols(A) == num_rows(B));
        C.resize(A.num_rows(), B.num_cols());
	ambient::push(ambient::gemm_diagonal_lhs_l,ambient::gemm_diagonal_lhs_c, A.get_data(), B, C);
        ambient::playout();
    }
    
    template<typename T>
    void svd(const p_dense_matrix<T>& A,
                   p_dense_matrix<T>& U,
                   p_dense_matrix<T>& V,
             typename associated_diagonal_matrix<p_dense_matrix<T> >::type& S)
    {
        BOOST_CONCEPT_ASSERT((blas::Matrix<p_dense_matrix<T> >));
	typename p_dense_matrix<T>::size_type k = std::min(num_rows(A), num_cols(A));
        U.resize(num_rows(A), k);
        V.resize(k, num_cols(A));
        S.resize(k, k);
	ambient::push(ambient::svd_l_scalapack, ambient::svd_c_scalapack, A, U, V, S.get_data());
        ambient::playout();
    }
    
    template<typename T>
    void qr(p_dense_matrix<T> M,
            p_dense_matrix<T> & Q,
            p_dense_matrix<T> & R)
    {
        assert(false);
        /* implement thin QR decomposition, i.e. for a (m,n) matrix, where m >= n, the result should be
         Q: (m,n)
         R: (n,n) */
    }

    template<typename T>
    void syev(p_dense_matrix<T> M,
              p_dense_matrix<T> & evecs,
              std::vector<double> & evals)
    {
        assert(false);
        assert(num_rows(M) == num_cols(M));
        assert(evals.size() == num_rows(M));
    //    boost::numeric::bindings::lapack::syevd('V', M, evals);    // <- push kernel
        // to be consistent with the SVD, I reorder in decreasing order
        std::reverse(evals.begin(), evals.end());
        // and the same with the matrix
        evecs.resize(num_rows(M), num_cols(M));
        for (std::size_t c = 0; c < num_cols(M); ++c)
			std::copy(column(M, c).first, column(M, c).second,
                      column(evecs, num_cols(M)-1-c).first);
        ambient::playout();
    }

    template<typename T>
    void syev(p_dense_matrix<T> M,
              p_dense_matrix<T> & evecs,
              typename associated_diagonal_matrix<p_dense_matrix<T> >::type & evals)
    {
        assert(false);

        assert(num_rows(M) == num_cols(M));
        std::vector<double> evals_(num_rows(M));
        syev(M, evecs, evals_);  
  //      evals = typename associated_diagonal_matrix<p_dense_matrix<T>::type(evals_); to modify
        ambient::playout();
    }
 
   template<typename T>
   void validation(const p_dense_matrix<T> & A_ambient, const p_dense_matrix<T> & B_scala)
   {
       ambient::push(ambient::validation_l, ambient::validation_c, A_ambient, B_scala);
        ambient::playout();
   }

   template<typename T, class Generator>
   void generate(p_dense_matrix<T>& A, Generator g)
   { // empty - just for compatibility
       A.set_init(ambient::random_i<T>);
   }

   template<class Matrix>
   void copy(typename associated_diagonal_matrix<Matrix>::type & S,typename associated_vector<Matrix>::type & allS)
   {
       assert(false);
   }
   
   template<class T>
   void sort(typename associated_vector<T>::type & allS)
   {
       assert(false);
   }

   template<class T>
   void reverse(typename associated_vector<T>::type & allS)
   {
       assert(false);
   }
 
} /* namespace blas */

#endif
