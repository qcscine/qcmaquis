#ifndef __ALPS_MATRIX_INTERFACE_HPP__
#define __ALPS_MATRIX_INTERFACE_HPP__

#include "matrix_concept_check.hpp"
#include "general_matrix.hpp"

namespace blas
{

// This macro creates free functions that call member functions with the same
// name, e.g. swap_columns(A,i,j) -> A.swap_columns(i,j)
#define COMMA ,
#define IMPLEMENT_FORWARDING(TEMPLATE_PARS,TYPE,RET,NAME,ARGS,VARS) \
template TEMPLATE_PARS \
RET NAME ARGS \
{ \
    /* BOOST_CONCEPT_ASSERT((blas::Matrix<TYPE>)); */ \
    return m.NAME VARS; \
} 

// num_rows(), num_columns(), swap_rows(), swap_columns()
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     typename general_matrix<T COMMA MemoryBlock>::size_type, num_rows, (general_matrix<T, MemoryBlock> const& m), () )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     typename general_matrix<T COMMA MemoryBlock>::size_type, num_columns, (general_matrix<T, MemoryBlock> const& m), () )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     void, swap_rows, (general_matrix<T, MemoryBlock>& m, typename general_matrix<T, MemoryBlock>::size_type i1, typename general_matrix<T, MemoryBlock>::size_type i2), (i1,i2) )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     void, swap_columns, (general_matrix<T, MemoryBlock>& m, typename general_matrix<T, MemoryBlock>::size_type i1, typename general_matrix<T, MemoryBlock>::size_type i2), (i1,i2) )
    
//
// Matrix Iterator Interface
//

#define ITERATOR_PAIR(TYPE, ITERATOR) \
std::pair<typename TYPE::ITERATOR, typename TYPE::ITERATOR>

IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(general_matrix<T COMMA MemoryBlock>, row_element_iterator), row,
                     (general_matrix<T COMMA MemoryBlock> & m,
                      typename general_matrix<T COMMA MemoryBlock>::size_type i),
                     (i) )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(general_matrix<T COMMA MemoryBlock>, const_row_element_iterator), row,
                     (general_matrix<T COMMA MemoryBlock> const& m,
                      typename general_matrix<T COMMA MemoryBlock>::size_type i),
                     (i) )    

IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(general_matrix<T COMMA MemoryBlock>, column_element_iterator), column,
                     (general_matrix<T COMMA MemoryBlock> & m,
                      typename general_matrix<T COMMA MemoryBlock>::size_type i),
                     (i) )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(general_matrix<T COMMA MemoryBlock>, const_column_element_iterator), column,
                     (general_matrix<T COMMA MemoryBlock> const& m,
                      typename general_matrix<T COMMA MemoryBlock>::size_type i),
                     (i) )  

IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(general_matrix<T COMMA MemoryBlock>, element_iterator), elements,
                     (general_matrix<T COMMA MemoryBlock>& m), () )

IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     ITERATOR_PAIR(general_matrix<T COMMA MemoryBlock>, const_element_iterator), elements,
                     (general_matrix<T COMMA MemoryBlock> const& m), () )
#undef ITERATOR_PAIR
#undef IMPLEMENT_FORWARDING
#undef COMMA
} //namespace blas

#endif //__ALPS_MATRIX_INTERFACE_HPP__
