#ifndef __ALPS_P_DIAGONAL_MATRIX_H
#define __ALPS_P_DIAGONAL_MATRIX_H


namespace blas {

    template<class FullMatrixClass>
    struct associated_p_diagonal_matrix { };
    
/**
idea the container is a p_dense_matrix of one row and n columns,
we need it to avoid dependency inside SVD kernels and other 
*/
    template<typename T>
    class p_diagonal_matrix
    {
    public:
        typedef T                       value_type;
        typedef T&                      reference;
        typedef T const&                const_reference;
        typedef size_t             size_type;
        typedef std::ptrdiff_t          difference_type;
        
        typedef typename std::vector<T>::iterator element_iterator;
        typedef typename std::vector<T>::const_iterator const_element_iterator;

	//to do, usefull ? check with SVD ....
       // template<class Vector>
      //  p_diagonal_matrix(Vector const & init)
      //  : data_(init.begin(), init.end()) { }
        p_diagonal_matrix(size_t rows = 0,  T const & init = T());
        size_t num_rows() const;
        size_t num_cols() const;
        T const & operator[](size_t i) const;
        T & operator[](size_t i); 
        T const & operator()(size_t i, size_t j) const;
        T & operator()(size_t i, size_t j);
        std::pair<element_iterator, element_iterator> elements();
        std::pair<const_element_iterator, const_element_iterator> elements() const;
        void remove_rows(size_t k, size_t n = 1);
        void remove_cols(size_t k, size_t n = 1);
        void resize(size_t rows, size_t cols, T v = T());
	template< class T1> 
        friend std::ostream & operator <<(std::ostream& os, p_diagonal_matrix<T1> const &m);
        const ambient::p_dense_matrix<T> & get_data() const; 
        ambient::p_dense_matrix<T> & get_data(); 
    private:
        ambient::p_dense_matrix<T> data_;
    };
}

#include "p_dense_matrix/p_diagonal_matrix.hpp"
#endif
