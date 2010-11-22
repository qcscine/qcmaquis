#ifndef ALPS_RESIZABLE_MATRIX_INTERFACE_HPP
#define ALPS_RESIZABLE_MATRIX_INTERFACE_HPP

#include "resizable_matrix_concept_check.hpp"

namespace blas
{

// This macro creates free functions that call member functions with the same
// name, e.g. swap_columns(A,i,j) -> A.swap_columns(i,j)
#define IMPLEMENT_FORWARDING(RET,NAME,ARGS,VARS) \
template <typename ResizableMatrix> \
RET NAME ARGS \
{ \
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<ResizableMatrix>)); \
    return m.NAME VARS; \
}

// resize(), remove_row(), remove_column()
IMPLEMENT_FORWARDING(void, resize, (ResizableMatrix& m, typename ResizableMatrix::size_type i, typename ResizableMatrix::size_type j), (i,j) )
IMPLEMENT_FORWARDING(void, resize, (ResizableMatrix& m, typename ResizableMatrix::size_type i, typename ResizableMatrix::size_type j, typename ResizableMatrix::value_type t), (i,j,t) )
IMPLEMENT_FORWARDING(void, remove_rows, (ResizableMatrix& m, typename ResizableMatrix::size_type i, typename ResizableMatrix::difference_type k = 1), (i,k) )
IMPLEMENT_FORWARDING(void, remove_columns, (ResizableMatrix& m, typename ResizableMatrix::size_type j, typename ResizableMatrix::difference_type k = 1), (j,k) )

#undef IMPLEMENT_FORWARDING


#define IMPLEMENT_ITER_FCT_FORWARDING(RET,NAME,ARGS,VARS) \
template<typename ResizableMatrix, typename InputIterator> \
RET NAME ARGS \
{ \
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<ResizableMatrix>)); \
    return m.NAME VARS; \
}

#define INPUT_ITERATOR_PAIR std::pair<InputIterator,InputIterator>

//append_row(), append_column(), insert_row(), insert_column()
IMPLEMENT_ITER_FCT_FORWARDING(void, append_rows, (ResizableMatrix& m, INPUT_ITERATOR_PAIR range, typename ResizableMatrix::difference_type k = 1), (range,k) )
IMPLEMENT_ITER_FCT_FORWARDING(void, append_columns, (ResizableMatrix& m, INPUT_ITERATOR_PAIR range, typename ResizableMatrix::difference_type k = 1), (range,k) )
IMPLEMENT_ITER_FCT_FORWARDING(void, insert_rows, (ResizableMatrix& m, typename ResizableMatrix::size_type i, INPUT_ITERATOR_PAIR range, typename ResizableMatrix::difference_type k = 1), (i,range,k) )
IMPLEMENT_ITER_FCT_FORWARDING(void, insert_columns, (ResizableMatrix& m, typename ResizableMatrix::size_type j, INPUT_ITERATOR_PAIR range, typename ResizableMatrix::difference_type k = 1), (j,range,k) )

#undef INPUT_ITERATOR_PAIR
#undef IMPLEMENT_ITER_FCT_FORWARDING

} // namespace blas

#endif //ALPS_RESIZABLE_MATRIX_INTERFACE_HPP
