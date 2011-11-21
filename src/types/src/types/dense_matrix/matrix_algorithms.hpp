#ifndef __ALPS_MATRIX_ALGORITHMS_HPP__
#define __ALPS_MATRIX_ALGORITHMS_HPP__

#include <types/dense_matrix/matrix_concept_check.hpp>

#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include <types/dense_matrix/dense_matrix.h>
#include <types/dense_matrix/diagonal_matrix.h> // remove?

#include <utils/timings.h>


namespace maquis {
    namespace types {   

    template <typename Matrix>
    Matrix transpose(Matrix const& m) 
    {
        static Timer timer("transpose");
        timer.begin();
        
        BOOST_CONCEPT_ASSERT((maquis::types::Matrix<Matrix>)); 
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
        BOOST_CONCEPT_ASSERT((maquis::types::Matrix<Matrix>)); 
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
 
    template<class Matrix>
    Matrix conjugate(Matrix M)
    {
        M.inplace_conjugate();
        return M;
    }

    template<class Matrix> Matrix join(Matrix const & a, Matrix const & b)
    {
        Matrix ret(num_rows(a)+num_rows(b), num_cols(a)+num_cols(b));
        
        typedef typename Matrix::size_type st;
        
        for (st r = 0; r < num_rows(a); ++r)
            for (st c = 0; c < num_cols(a); ++c)
                ret(r, c) = a(r, c);
        
        for (st r = 0; r < num_rows(b); ++r)
            for (st c = 0; c < num_cols(b); ++c)
                ret(r+num_rows(a), c+num_cols(a)) = b(r, c);
        
        return ret;
    }

    } // namespace types
} //namespace maquis

#endif //__ALPS_MATRIX_ALGORITHMS_HPP__
