#ifndef __ALPS_P_MATRIX_INTERFACE_HPP__
#define __ALPS_P_MATRIX_INTERFACE_HPP__

#include "p_dense_matrix/concept/matrix_concept_check.hpp"
#include "p_dense_matrix/p_dense_matrix.h"

namespace blas
{

// This macro creates free functions that call member functions with the same name
#define COMMA ,
#define IMPLEMENT_FORWARDING(TEMPLATE_PARS,TYPE,RET,NAME,ARGS,VARS) \
template TEMPLATE_PARS \
RET NAME ARGS \
{ \
    return m.NAME VARS; \
} 


//    BOOST_CONCEPT_ASSERT((blas::Matrix<TYPE>)); this line was before return, but pb with boost 


// num_rows(), num_cols()
IMPLEMENT_FORWARDING(<typename T>, p_dense_matrix<T>,
                     typename p_dense_matrix<T>::size_type, num_rows, (p_dense_matrix<T> const& m), () )
IMPLEMENT_FORWARDING(<typename T>, p_dense_matrix<T>,
                     typename p_dense_matrix<T>::size_type, num_cols, (p_dense_matrix<T> const& m), () )
//
// Matrix Iterator Interface
//

#define ITERATOR_PAIR(TYPE, ITERATOR) \
std::pair<typename TYPE::ITERATOR, typename TYPE::ITERATOR>

IMPLEMENT_FORWARDING(<typename T>, p_dense_matrix<T>,
                     ITERATOR_PAIR(p_dense_matrix<T>, row_element_iterator), row,
                     (p_dense_matrix<T> & m,
                      typename p_dense_matrix<T>::size_type i),
                     (i) )
IMPLEMENT_FORWARDING(<typename T>, p_dense_matrix<T>,
                     ITERATOR_PAIR(p_dense_matrix<T>, const_row_element_iterator), row,
                     (p_dense_matrix<T> const& m,
                      typename p_dense_matrix<T>::size_type i),
                     (i) )    

IMPLEMENT_FORWARDING(<typename T>, p_dense_matrix<T>,
                     ITERATOR_PAIR(p_dense_matrix<T>, column_element_iterator), column,
                     (p_dense_matrix<T> & m,
                      typename p_dense_matrix<T>::size_type i),
                     (i) )
IMPLEMENT_FORWARDING(<typename T>, p_dense_matrix<T>,
                     ITERATOR_PAIR(p_dense_matrix<T>, const_column_element_iterator), column,
                     (p_dense_matrix<T> const& m,
                      typename p_dense_matrix<T>::size_type i),
                     (i) )  

IMPLEMENT_FORWARDING(<typename T>, p_dense_matrix<T>,
                     ITERATOR_PAIR(p_dense_matrix<T>, element_iterator), elements,
                     (p_dense_matrix<T>& m), () )

IMPLEMENT_FORWARDING(<typename T>, p_dense_matrix<T>,
                     ITERATOR_PAIR(p_dense_matrix<T>, const_element_iterator), elements,
                     (p_dense_matrix<T> const& m), () )
#undef ITERATOR_PAIR
#undef IMPLEMENT_FORWARDING
#undef COMMA
} //namespace blas

#endif //__ALPS_MATRIX_INTERFACE_HPP__
