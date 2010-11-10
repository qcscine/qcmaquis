#ifndef __ALPS_MATRIX_INTERFACE_HPP__
#define __ALPS_MATRIX_INTERFACE_HPP__

#include "matrix_concept_check.hpp"

namespace blas
{

#define IMPLEMENT_FORWARDING(RET,NAME,ARGS,VARS) \
template <typename MatrixType> \
RET NAME ARGS \
{ \
    BOOST_CONCEPT_ASSERT((blas::Matrix<MatrixType>)); \
    return m.NAME VARS; \
} \

// num_rows(), num_columns(), swap_rows(), swap_columns()
IMPLEMENT_FORWARDING(typename MatrixType::size_type, num_rows, (MatrixType const& m), () )
IMPLEMENT_FORWARDING(typename MatrixType::size_type, num_columns, (MatrixType const& m), () )
IMPLEMENT_FORWARDING(void, swap_rows, (MatrixType& m, typename MatrixType::size_type i1, typename MatrixType::size_type i2), (i1,i2) )
IMPLEMENT_FORWARDING(void, swap_columns, (MatrixType& m, typename MatrixType::size_type i1, typename MatrixType::size_type i2), (i1,i2) )

//
// Matrix Iterator Interface
//

#define ITERATOR_PAIR(ITERATOR) \
std::pair<typename MatrixType::ITERATOR, typename MatrixType::ITERATOR>

IMPLEMENT_FORWARDING( ITERATOR_PAIR(row_element_iterator), row, (MatrixType& m, typename MatrixType::size_type i), (i) )
IMPLEMENT_FORWARDING( ITERATOR_PAIR(const_row_element_iterator), row, (MatrixType const& m, typename MatrixType::size_type i), (i) )
IMPLEMENT_FORWARDING( ITERATOR_PAIR(column_element_iterator), column, (MatrixType& m, typename MatrixType::size_type j), (j) )
IMPLEMENT_FORWARDING( ITERATOR_PAIR(const_column_element_iterator), column, (MatrixType const& m, typename MatrixType::size_type j), (j) )
#undef ITERATOR_PAIR
#undef IMPLEMENT_FORWARDING
} //namespace blas

#endif //__ALPS_MATRIX_INTERFACE_HPP__
