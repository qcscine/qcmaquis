#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

#include "p_dense_matrix/p_diagonal_matrix.h"


namespace blas {
  
//to do wrapper scala first
    template<typename T, class Matrix>
    void gemm(Matrix const & m1, diagonal_matrix<T> const & m2, Matrix & m3)
    {

	// to do scalapack kernel pdgemv and home made kernel 
        assert(num_columns(m1) == num_rows(m2));
        resize(m3, num_rows(m1), num_columns(m2));
        for (std::size_t i = 0; i < num_rows(m1); ++i)
            for (std::size_t j = 0; j < num_columns(m2); ++j)
                m3(i,j) = m1(i,j) * m2(j,j);
    }
    
//to do wrapper scala first
    template<typename T, class Matrix>
    void gemm(diagonal_matrix<T> const & m1, Matrix const & m2, Matrix & m3)
    {
	// to do scalapack kernel pdgemv and home made kernel 
        assert(num_columns(m1) == num_rows(m2));
        resize(m3, num_rows(m1), num_columns(m2));
        for (std::size_t i = 0; i < num_rows(m1); ++i)
            for (std::size_t j = 0; j < num_columns(m2); ++j)
                m3(i,j) = m1(i,i) * m2(i,j);
    }

    template<typename T>
    typename diagonal_matrix<T>::size_type num_rows(diagonal_matrix<T> const & m)
    {
        return m.num_rows();
    }
    
    template<typename T>
    typename diagonal_matrix<T>::size_type num_columns(diagonal_matrix<T> const & m)
    {
        return m.num_columns();
    }
 

//to do pinned kernel
    template<typename T>
    diagonal_matrix<T> sqrt(diagonal_matrix<T> m)
    {
        std::transform(m.elements().first, m.elements().second, m.elements().first, utils::functor_sqrt());
        return m;
    }
   
/* become useless 
    template<typename T>
    std::ostream& operator<<(std::ostream& os, diagonal_matrix<T> const & m)
    {
        std::copy(m.elements().first, m.elements().second, std::ostream_iterator<T>(os, " "));
        return os;
    }
 */ 
 
    template<typename T>
    void remove_rows(diagonal_matrix<T> & m, std::size_t k, std::size_t n = 1)
    {
        m.remove_rows(k, n);
    }
    
    template<typename T>
    void remove_columns(diagonal_matrix<T> & m, std::size_t k, std::size_t n = 1)
    {
        m.remove_columns(k, n);
    }
   

//to check
    template<typename T>
    void resize(diagonal_matrix<T> & m, std::size_t r, std::size_t c, T v = T())
    {
        m.resize(r, c, v);
    }
    
    template<typename T>
    std::pair<typename diagonal_matrix<T>::element_iterator, typename diagonal_matrix<T>::element_iterator>
    elements(diagonal_matrix<T> & m)
    {
        return m.elements();
    }

    template<typename T>
    std::pair<typename diagonal_matrix<T>::const_element_iterator, typename diagonal_matrix<T>::const_element_iterator>
    elements(diagonal_matrix<T> const & m)
    {
        return m.elements();
    }
}

#endif
