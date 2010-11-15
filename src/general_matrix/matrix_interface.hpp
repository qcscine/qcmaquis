#ifndef __ALPS_MATRIX_INTERFACE_HPP__
#define __ALPS_MATRIX_INTERFACE_HPP__

#include "matrix_concept_check.hpp"

namespace blas
{

// This macro creates free functions that call member functions with the same
// name, e.g. swap_columns(A,i,j) -> A.swap_columns(i,j)
#define IMPLEMENT_FORWARDING(RET,NAME,ARGS,VARS) \
template <typename Matrix> \
RET NAME ARGS \
{ \
    BOOST_CONCEPT_ASSERT((blas::Matrix<Matrix>)); \
    return m.NAME VARS; \
} 

// num_rows(), num_columns(), swap_rows(), swap_columns()
IMPLEMENT_FORWARDING(typename Matrix::size_type, num_rows, (Matrix const& m), () )
IMPLEMENT_FORWARDING(typename Matrix::size_type, num_columns, (Matrix const& m), () )
IMPLEMENT_FORWARDING(void, swap_rows, (Matrix& m, typename Matrix::size_type i1, typename Matrix::size_type i2), (i1,i2) )
IMPLEMENT_FORWARDING(void, swap_columns, (Matrix& m, typename Matrix::size_type i1, typename Matrix::size_type i2), (i1,i2) )

//
// Matrix Iterator Interface
//

#define ITERATOR_PAIR(ITERATOR) \
std::pair<typename Matrix::ITERATOR, typename Matrix::ITERATOR>

IMPLEMENT_FORWARDING( ITERATOR_PAIR(row_element_iterator), row, (Matrix& m, typename Matrix::size_type i), (i) )
IMPLEMENT_FORWARDING( ITERATOR_PAIR(const_row_element_iterator), row, (Matrix const& m, typename Matrix::size_type i), (i) )
IMPLEMENT_FORWARDING( ITERATOR_PAIR(column_element_iterator), column, (Matrix& m, typename Matrix::size_type j), (j) )
IMPLEMENT_FORWARDING( ITERATOR_PAIR(const_column_element_iterator), column, (Matrix const& m, typename Matrix::size_type j), (j) )
#undef ITERATOR_PAIR
#undef IMPLEMENT_FORWARDING
} //namespace blas

#endif //__ALPS_MATRIX_INTERFACE_HPP__
