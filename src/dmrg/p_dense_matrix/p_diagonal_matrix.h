#ifndef __ALPS_P_DIAGONAL_MATRIX_H
#define __ALPS_P_DIAGONAL_MATRIX_H


namespace blas {

    template<class FullMatrixClass>
    struct associated_diagonal_matrix { };
    
/**
idea the container is a p_dense_matrix of one row and n columns,
we need it to avoid dependency inside SVD kernels and other 
*/
    template<typename T>
    class p_diagonal_matrix
    {
    public:
        typedef typename ambient::p_dense_matrix<T>          container; // can be ::dynamic
        typedef T                                            value_type;
        typedef T&                                           reference;
        typedef T const&                                     const_reference;
        typedef size_t                                       size_type;
        typedef size_t                                       difference_type;
        
        typedef typename container::element_iterator element_iterator;
        typedef typename container::const_element_iterator const_element_iterator;

      
	//to do, usefull ? check with SVD, yes due to the structure of pdgsvd 
       // template<class Vector>
      //  p_diagonal_matrix(Vector const & init)
      //  : data_(init.begin(), init.end()) { }
        p_diagonal_matrix(size_t rows = 0,  T const & init = T());
      //  p_diagonal_matrix(size_t rows,  T* const &  array);
        size_t num_rows() const;
        size_t num_cols() const;
        p_diagonal_matrix<T>& operator = (const p_diagonal_matrix<T>& rhs);
        T const & operator[](size_t i) const;
        T & operator[](size_t i); 
        T const & operator()(size_t i, size_t j) const;
        T & operator()(size_t i, size_t j);
        std::pair<element_iterator, element_iterator> elements();
        std::pair<const_element_iterator, const_element_iterator> elements() const;
        void remove_rows(size_t k, size_t n = 1);
        void remove_cols(size_t k, size_t n = 1);
        void resize(size_t rows, size_t cols, T v = T());
	template< class T1 > 
        friend std::ostream & operator <<(std::ostream& os, p_diagonal_matrix<T1> const &m);
        const container & get_data() const; 
        container & get_data();    
        element_iterator begin();
        const_element_iterator begin() const;
        element_iterator end();
        const_element_iterator end() const;
        std::size_t size() const;
/*
        void push_back(element_iterator);
        void push_back(element_iterator) const;
 */
   private:
        container data_;
    };
}

#include "p_dense_matrix/p_diagonal_matrix.hpp"
#endif
