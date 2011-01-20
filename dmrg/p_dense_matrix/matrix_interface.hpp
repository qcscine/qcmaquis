#ifndef __ALPS_MATRIX_INTERFACE_HPP__
#define __ALPS_MATRIX_INTERFACE_HPP__

#include "p_dense_matrix/matrix_concept_check.hpp"
#include "p_dense_matrix/p_dense_matrix.h"

namespace blas
{

// This macro creates free functions that call member functions with the same
// name, e.g. swap_columns(A,i,j) -> A.swap_columns(i,j)
#define COMMA ,
#define IMPLEMENT_FORWARDING(TEMPLATE_PARS,TYPE,RET,NAME,ARGS,VARS) \
template TEMPLATE_PARS \
RET NAME ARGS \
{ \
    BOOST_CONCEPT_ASSERT((blas::Matrix<TYPE>)); \
    return m.NAME VARS; \
} 

// num_rows(), num_columns(), swap_rows(), swap_columns()
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     typename p_dense_matrix<T COMMA MemoryBlock>::size_type, num_rows, (p_dense_matrix<T, MemoryBlock> const& m), () )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     typename p_dense_matrix<T COMMA MemoryBlock>::size_type, num_columns, (p_dense_matrix<T, MemoryBlock> const& m), () )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     void, swap_rows, (p_dense_matrix<T, MemoryBlock>& m, typename p_dense_matrix<T, MemoryBlock>::size_type i1, typename p_dense_matrix<T, MemoryBlock>::size_type i2), (i1,i2) )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     void, swap_columns, (p_dense_matrix<T, MemoryBlock>& m, typename p_dense_matrix<T, MemoryBlock>::size_type i1, typename p_dense_matrix<T, MemoryBlock>::size_type i2), (i1,i2) )
    
//
// Matrix Iterator Interface
//

#define ITERATOR_PAIR(TYPE, ITERATOR) \
std::pair<typename TYPE::ITERATOR, typename TYPE::ITERATOR>

IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(p_dense_matrix<T COMMA MemoryBlock>, row_element_iterator), row,
                     (p_dense_matrix<T COMMA MemoryBlock> & m,
                      typename p_dense_matrix<T COMMA MemoryBlock>::size_type i),
                     (i) )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(p_dense_matrix<T COMMA MemoryBlock>, const_row_element_iterator), row,
                     (p_dense_matrix<T COMMA MemoryBlock> const& m,
                      typename p_dense_matrix<T COMMA MemoryBlock>::size_type i),
                     (i) )    

IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(p_dense_matrix<T COMMA MemoryBlock>, column_element_iterator), column,
                     (p_dense_matrix<T COMMA MemoryBlock> & m,
                      typename p_dense_matrix<T COMMA MemoryBlock>::size_type i),
                     (i) )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(p_dense_matrix<T COMMA MemoryBlock>, const_column_element_iterator), column,
                     (p_dense_matrix<T COMMA MemoryBlock> const& m,
                      typename p_dense_matrix<T COMMA MemoryBlock>::size_type i),
                     (i) )  

IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(p_dense_matrix<T COMMA MemoryBlock>, element_iterator), elements,
                     (p_dense_matrix<T COMMA MemoryBlock>& m), () )

IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, p_dense_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(p_dense_matrix<T COMMA MemoryBlock>, const_element_iterator), elements,
                     (p_dense_matrix<T COMMA MemoryBlock> const& m), () )
#undef ITERATOR_PAIR
#undef IMPLEMENT_FORWARDING
#undef COMMA
} //namespace blas

#endif //__ALPS_MATRIX_INTERFACE_HPP__
