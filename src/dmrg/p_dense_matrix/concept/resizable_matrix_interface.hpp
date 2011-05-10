#ifndef __ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__
#define __ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__

#include "p_dense_matrix/concept/resizable_matrix_concept_check.hpp"

namespace blas
{

// resize(), remove_row(), remove_column()
template <typename T>
void resize(p_dense_matrix<T>& m, typename p_dense_matrix<T>::size_type i, typename p_dense_matrix<T>::size_type j)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T> >));
    return m.resize(i,j);
}

template <typename T>
void resize( p_dense_matrix<T>& m,
        typename p_dense_matrix<T>::size_type i,
        typename p_dense_matrix<T>::size_type j,
        typename p_dense_matrix<T>::value_type const& t )
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T> >));
    return m.resize(i,j); // maybe will add t back
}

template <typename T>
void remove_rows( p_dense_matrix<T>& m,
        typename p_dense_matrix<T>::size_type i,
        typename p_dense_matrix<T>::difference_type k = 1)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T> >));
    return m.remove_rows(i,k);
}

template <typename T>
void remove_cols( p_dense_matrix<T>& m,
        typename p_dense_matrix<T>::size_type j,
        typename p_dense_matrix<T>::difference_type k = 1)
{
    BOOST_CONCEPT_ASSERT((blas::ResizableMatrix<p_dense_matrix<T> >));
    return m.remove_cols(j,k);
}

} // namespace blas

#endif //__ALPS_RESIZABLE_MATRIX_INTERFACE_HPP__
