#ifndef __ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__
#define __ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__

#include "p_dense_matrix/resizable_matrix_concept_check.hpp"

namespace blas
{

// resize(), remove_row(), remove_column()
template <typename T, typename MemoryBlock>
void resize(p_dense_matrix<T,MemoryBlock>& m, typename p_dense_matrix<T,MemoryBlock>::size_type i, typename p_dense_matrix<T,MemoryBlock>::size_type j)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T,MemoryBlock> >));
    return m.resize(i,j);
}

template <typename T, typename MemoryBlock>
void resize( p_dense_matrix<T,MemoryBlock>& m,
        typename p_dense_matrix<T,MemoryBlock>::size_type i,
        typename p_dense_matrix<T,MemoryBlock>::size_type j,
        typename p_dense_matrix<T,MemoryBlock>::value_type const& t )
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T,MemoryBlock> >));
    return m.resize(i,j,t);
}

template <typename T, typename MemoryBlock>
void remove_rows( p_dense_matrix<T,MemoryBlock>& m,
        typename p_dense_matrix<T,MemoryBlock>::size_type i,
        typename p_dense_matrix<T,MemoryBlock>::difference_type k = 1)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T,MemoryBlock> >));
    return m.remove_rows(i,k);
}

template <typename T, typename MemoryBlock>
void remove_columns( p_dense_matrix<T,MemoryBlock>& m,
        typename p_dense_matrix<T,MemoryBlock>::size_type j,
        typename p_dense_matrix<T,MemoryBlock>::difference_type k = 1)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T,MemoryBlock> >));
    return m.remove_columns(j,k);
}

//append_row(), append_column(), insert_row(), insert_column()
#define INPUT_ITERATOR_PAIR std::pair<InputIterator,InputIterator>

template <typename T, typename MemoryBlock, typename InputIterator>
void append_rows( p_dense_matrix<T,MemoryBlock>& m, INPUT_ITERATOR_PAIR range,
        typename p_dense_matrix<T,MemoryBlock>::difference_type k = 1)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T,MemoryBlock> >));
    return m.append_rows(range,k);
}

template <typename T, typename MemoryBlock, typename InputIterator>
void append_columns( p_dense_matrix<T,MemoryBlock>& m, INPUT_ITERATOR_PAIR range,
        typename p_dense_matrix<T,MemoryBlock>::difference_type k = 1)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T,MemoryBlock> >));
    return m.append_columns(range,k);
}

template <typename T, typename MemoryBlock, typename InputIterator>
void insert_rows( p_dense_matrix<T,MemoryBlock>& m,
        typename p_dense_matrix<T,MemoryBlock>::size_type i,
        INPUT_ITERATOR_PAIR range,
        typename p_dense_matrix<T,MemoryBlock>::difference_type k = 1)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T,MemoryBlock> >));
    return m.insert_rows(i,range,k);
}

template <typename T, typename MemoryBlock, typename InputIterator>
void insert_columns( p_dense_matrix<T,MemoryBlock>& m,
        typename p_dense_matrix<T,MemoryBlock>::size_type j,
        INPUT_ITERATOR_PAIR range,
        typename p_dense_matrix<T,MemoryBlock>::difference_type k = 1)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T,MemoryBlock> >));
    return m.insert_columns(j,range,k);
}

#undef INPUT_ITERATOR_PAIR

} // namespace blas

#endif //__ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__
