#ifndef __ALPS_DENSE_MATRIX_ALGORITHMS_HPP__
#define __ALPS_DENSE_MATRIX_ALGORITHMS_HPP__

#include "p_dense_matrix/concept/matrix_concept_check.hpp"

#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include "p_dense_matrix/p_dense_matrix.h"
#include "utils/timings.h"

namespace blas
{
    namespace detail {
        template<typename T> struct sv_type { typedef T type; };
        template<typename T>
        struct sv_type<std::complex<T> > { typedef T type; };
    }

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
        printf("Attempting to perform INPCONJ\n");
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
        Matrix ret(size, size);
        for (typename Matrix::size_type k = 0; k < size; ++k)
            ret(k,k) = 1;
        return ret;
    }

    template<typename T>
    void gemm( const p_dense_matrix<T>  & A, const p_dense_matrix<T> & B, p_dense_matrix<T> & C)
    {
        printf("Attempting to perform GEMM\n");
	ambient::push(ambient::gemm_l_scalapack_kernel, ambient::gemm_c_scalapack_kernel, A,B,C );
    }
    
    template<typename T>
    void svd(const p_dense_matrix<T> &  M,
             p_dense_matrix<T> & U,
             p_dense_matrix<T>& V,
             typename associated_diagonal_matrix<p_dense_matrix<T> >::type & S)
    {
        //BOOST_CONCEPT_ASSERT((blas::Matrix<p_dense_matrix<T> >));
        printf("Attempting to perform SVD\n");
	typename p_dense_matrix<T>::size_type k = std::min(num_rows(M), num_cols(M));
        U.resize(num_rows(M), k);
        V.resize(k, num_cols(M));
       // std::vector<typename detail::sv_type<T>::type> S_(k);
	//const T * S_= new T[k];
	T* const *  S_ = new T*((T*)malloc(sizeof(T)*k));
	ambient::push(ambient::null_l_scalapack_svd_kernel, ambient::svd_c_scalapack_kernel, M, U, V, *S_);
        ambient::playout();
        p_diagonal_matrix<T> R(k,*S_) ; 
        S.get_data() = R.get_data(); //typename associated_p_diagonal_matrix<p_dense_matrix<T> >::type (R); 
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
    
} /* namespace blas */

#endif
