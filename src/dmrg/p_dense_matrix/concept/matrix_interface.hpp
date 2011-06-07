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

#undef IMPLEMENT_FORWARDING
#undef COMMA
} //namespace blas

#endif //__ALPS_MATRIX_INTERFACE_HPP__
