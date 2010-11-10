#ifndef __ALPS_MATRIX_INTERFACE_HPP__
#define __ALPS_MATRIX_INTERFACE_HPP__

#include "matrix_concept_check.hpp"
#include <iostream>

namespace blas
{

template <typename MatrixType>
typename MatrixType::size_type num_rows(MatrixType const& m)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    return m.num_rows();
}

template <typename MatrixType>
typename MatrixType::size_type num_columns(MatrixType const& m)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    return m.num_columns();
}

template <typename MatrixType>
void swap_rows(MatrixType& m, typename MatrixType::size_type i1, typename MatrixType::size_type i2)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    m.swap_rows(i1,i2);
}

template <typename MatrixType>
void swap_columns(MatrixType& m, typename MatrixType::size_type i1, typename MatrixType::size_type i2)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    m.swap_columns(i1,i2);
}


//
// Matrix Iterator Interface
//

template <typename MatrixType>
std::pair<typename MatrixType::row_element_iterator, typename MatrixType::row_element_iterator> row(MatrixType& m, typename MatrixType::size_type i)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    return m.row(i);
}

template <typename MatrixType>
std::pair<typename MatrixType::const_row_element_iterator, typename MatrixType::const_element_row_iterator> row(MatrixType const& m, typename MatrixType::size_type i)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    return m.row(i);
}

template <typename MatrixType>
std::pair<typename MatrixType::column_element_iterator, typename MatrixType::column_element_iterator> column(MatrixType& m, typename MatrixType::size_type j)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    return m.column(j);
}

template <typename MatrixType>
std::pair<typename MatrixType::const_column_element_iterator, typename MatrixType::const_column_element_iterator> column(MatrixType const& m, typename MatrixType::size_type j)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    return m.column(j);
}

//
// Operators
//


template <typename MatrixType>
boost::enable_if<
const MatrixType operator + (MatrixType a, MatrixType const& b)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    a += b;
    return a;
}

template <typename MatrixType>
const MatrixType operator - (MatrixType a, MatrixType const& b)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    a -= b;
    return a;
}

template <typename MatrixType>
const MatrixType operator * (MatrixType m, typename MatrixType::value_type const& t)
{
    return m*=t;
}

template <typename MatrixType>
const MatrixType operator * (typename MatrixType::value_type const& t, MatrixType m)
{
    return m*=t;
}
    
template<typename MatrixType, typename VectorType>
const VectorType operator * (MatrixType const& m, VectorType const& v)
{
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
    //TODO BOOST_CONCEPT_ASSERT((blas::Vector<VectorType>));
    assert( num_columns(m) == v.size() );
    VectorType result(num_rows(m));
    // Simple Matrix * Vector
    for(typename MatrixType::size_type i = 0; i < num_rows(m); ++i)
    {
        for(typename MatrixType::size_type j=0; j < num_columns(m); ++j)
        {
            result(i) = m(i,j) * v(j);
        }
    }
    return result;
}
    
template<typename MatrixType>
const MatrixType operator * (MatrixType const& m1, MatrixType const& m2)
{
    using detail::matrix_matrix_multiply;
    return matrix_matrix_multiply(m1,m2);
}

    
//    
//template <typename MatrixType>
//std::ostream& operator << (std::ostream& o, MatrixType const& m)
//{
//    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>));
//    for(typename MatrixType::size_type i=0; i< num_rows(m); ++i)
//    {
//        for(typename MatrixType::size_type j=0; j < num_columns(m); ++j)
//        {
//            typename MatrixType::value_type value(m(i,j));
//            o<<value<<" ";
//        }
//        o<<std::endl;
//    }
//    return o;
//}
//
} //namespace blas

#endif //__ALPS_MATRIX_INTERFACE_HPP__
