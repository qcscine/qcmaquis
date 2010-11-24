#ifndef __ALPS_MATRIX_ALGORITHMS_HPP__
#define __ALPS_MATRIX_ALGORITHMS_HPP__

#include "matrix_concept_check.hpp"

#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include "general_matrix.hpp"
#include "diagonal_matrix.h"

namespace blas
{
    template <typename Matrix>
    Matrix transpose(Matrix const& m) 
    {
        BOOST_CONCEPT_ASSERT((blas::Matrix<Matrix>)); 
        // TODO: perhaps this could return a proxy object
        Matrix tmp(num_columns(m), num_rows(m));
        for(typename Matrix::size_type i=0; i < num_rows(m); ++i){
            for(typename Matrix::size_type j=0; j < num_columns(m); ++j){
                tmp(j,i) = m(i,j);
            }
        }
        return tmp;
    }
    
    template <typename Matrix>
    const typename Matrix::value_type trace(Matrix const& m)
    {
        BOOST_CONCEPT_ASSERT((blas::Matrix<Matrix>)); 
        assert(num_rows(m) == num_columns(m));
        typename Matrix::value_type tr(m(0,0));
        for(typename Matrix::size_type i = 1; i<num_rows(m); ++i)
            tr += m(i,i);
        return tr;
    }

    template<typename T, class MemoryBlock>
    void svd(general_matrix<T, MemoryBlock> & M,
             general_matrix<T, MemoryBlock> & U,
             general_matrix<T, MemoryBlock>& V,
             diagonal_matrix<T> & S)
    {
        BOOST_CONCEPT_ASSERT((blas::Matrix<general_matrix<T, MemoryBlock> >));
        typename general_matrix<T, MemoryBlock>::size_type k = std::min(num_rows(M), num_columns(M));
        resize(U, num_rows(M), k);
        resize(V, k, num_columns(M));
        
        std::vector<double> S_(k);
        boost::numeric::bindings::lapack::gesdd('S', M, S_, U, V);
        
        S = diagonal_matrix<T>(S_);
    }
}

#endif //__ALPS_MATRIX_ALGORITHMS_HPP__
