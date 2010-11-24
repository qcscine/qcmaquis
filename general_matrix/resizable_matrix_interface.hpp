#ifndef __ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__
#define __ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__

#include "resizable_matrix_concept_check.hpp"

namespace blas
{

// This macro creates free functions that call member functions with the same
// name, e.g. swap_columns(A,i,j) -> A.swap_columns(i,j)
#define COMMA ,
#define IMPLEMENT_FORWARDING(TEMPLATE_PARS,TYPE,RET,NAME,ARGS,VARS) \
template TEMPLATE_PARS \
RET NAME ARGS \
{ \
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<TYPE>)); \
    return m.NAME VARS; \
}

// resize(), remove_row(), remove_column()
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     void, resize,
                     (general_matrix<T COMMA MemoryBlock>& m,
                      typename general_matrix<T COMMA MemoryBlock>::size_type i,
                      typename general_matrix<T COMMA MemoryBlock>::size_type j),
                     (i,j) )
IMPLEMENT_FORWARDING(<typename T COMMA class MemoryBlock>, general_matrix<T COMMA MemoryBlock>,
                     void, resize,
                     (general_matrix<T COMMA MemoryBlock>& m,
                      typename general_matrix<T COMMA MemoryBlock>::size_type i,
                      typename general_matrix<T COMMA MemoryBlock>::size_type j,
                      typename general_matrix<T COMMA MemoryBlock>::value_type t),
                     (i,j, t) )
    
#warning FIXME
//IMPLEMENT_FORWARDING(void, remove_rows, (ResizableMatrix& m, typename ResizableMatrix::size_type i, typename ResizableMatrix::difference_type k = 1), (i,k) )
//IMPLEMENT_FORWARDING(void, remove_columns, (ResizableMatrix& m, typename ResizableMatrix::size_type j, typename ResizableMatrix::difference_type k = 1), (j,k) )

#undef IMPLEMENT_FORWARDING


#define IMPLEMENT_ITER_FCT_FORWARDING(RET,NAME,ARGS,VARS) \
template<typename ResizableMatrix, typename InputIterator> \
RET NAME ARGS \
{ \
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<ResizableMatrix>)); \
    return m.NAME VARS; \
}

#define INPUT_ITERATOR_PAIR std::pair<InputIterator,InputIterator>

#warning FIXME
//append_row(), append_column(), insert_row(), insert_column()
IMPLEMENT_ITER_FCT_FORWARDING(void, append_rows, (ResizableMatrix& m, INPUT_ITERATOR_PAIR range, typename ResizableMatrix::difference_type k = 1), (range,k) )
IMPLEMENT_ITER_FCT_FORWARDING(void, append_columns, (ResizableMatrix& m, INPUT_ITERATOR_PAIR range, typename ResizableMatrix::difference_type k = 1), (range,k) )
IMPLEMENT_ITER_FCT_FORWARDING(void, insert_rows, (ResizableMatrix& m, typename ResizableMatrix::size_type i, INPUT_ITERATOR_PAIR range, typename ResizableMatrix::difference_type k = 1), (i,range,k) )
IMPLEMENT_ITER_FCT_FORWARDING(void, insert_columns, (ResizableMatrix& m, typename ResizableMatrix::size_type j, INPUT_ITERATOR_PAIR range, typename ResizableMatrix::difference_type k = 1), (j,range,k) )

#undef INPUT_ITERATOR_PAIR
#undef IMPLEMENT_ITER_FCT_FORWARDING
#undef COMMA

} // namespace blas

#endif //__ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__
