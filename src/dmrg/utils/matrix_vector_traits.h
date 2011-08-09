
#ifndef __ALPS_MATRIX_VECTOR_TRAITS_HPP__
#define __ALPS_MATRIX_VECTOR_TRAITS_HPP__

#include <complex>

#include "utils/utils.hpp"

namespace blas {
    
    class Transpose { };
    class NoTranspose { };
    
    template<class M>
    struct associated_diagonal_matrix { };
    
    template<class M>
    struct associated_vector { };
    
    template<class M>
    struct associated_real_diagonal_matrix { };
    
    template<class M>
    struct associated_real_vector { };
    
    namespace detail {
        template<class T> struct real_type { typedef T type; };
        template<class T> struct real_type<std::complex<T> > { typedef T type; };
    }
    
    namespace detail
    {
        template<class Tag> struct evaluate_tag;
        template<> struct evaluate_tag<blas::NoTranspose>
        {
            template<class Matrix>
            static Matrix const & eval(Matrix const & m) { return m; }
        };
        template<> struct evaluate_tag<blas::Transpose>
        {
            template<class Matrix>
            static Matrix eval(Matrix const & m) { return transpose(m); }
        };
        
        template<class Tag> struct dims;
        
        template<> struct dims<NoTranspose>
        {
            template<class Matrix>
            static typename Matrix::size_type
            first(Matrix const & m)
            { return num_rows(m); }
            
            template<class Matrix>
            static typename Matrix::size_type
            second(Matrix const & m)
            { return num_cols(m); }
        };
        
        template<> struct dims<Transpose>
        {
            template<class Matrix>
            static typename Matrix::size_type
            first(Matrix const & m)
            { return num_cols(m); }
            
            template<class Matrix>
            static typename Matrix::size_type
            second(Matrix const & m)
            { return num_rows(m); }
        };
    }
}

template<class Matrix1, class Matrix2, class Matrix3, class Tag1, class Tag2>
void gemm(Matrix1 const & A, Tag1,
          Matrix2 const & B, Tag2,
          Matrix3 & C)
{
    gemm(blas::detail::evaluate_tag<Tag1>::eval(A),
         blas::detail::evaluate_tag<Tag2>::eval(B),
         C);
}

template<class Matrix1, class Matrix2, class Tag1, class Tag2>
std::pair<std::size_t, std::size_t> result_size(Matrix1 const & A, Tag1,
                                                Matrix2 const & B, Tag2)
{
    return std::make_pair(blas::detail::dims<Tag1>::first(A), blas::detail::dims<Tag2>::second(B));
}

#endif
