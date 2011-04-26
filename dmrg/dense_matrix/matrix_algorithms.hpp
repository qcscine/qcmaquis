#ifndef __ALPS_MATRIX_ALGORITHMS_HPP__
#define __ALPS_MATRIX_ALGORITHMS_HPP__

#include "dense_matrix/matrix_concept_check.hpp"

#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/diagonal_matrix.h" // remove?

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
        return Matrix::identity_matrix(size);
    }
}

#endif //__ALPS_MATRIX_ALGORITHMS_HPP__
