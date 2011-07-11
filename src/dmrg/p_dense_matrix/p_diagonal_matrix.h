#ifndef __ALPS_P_DIAGONAL_MATRIX_H
#define __ALPS_P_DIAGONAL_MATRIX_H

namespace blas {
    
    template<typename T>
    class p_diagonal_matrix
    {
    public:
        typedef typename ambient::p_dense_matrix<T> container; // can be ::dynamic
        typedef T                                   value_type;
        typedef T&                                  reference;
        typedef const T&                            const_reference;
        typedef size_t                              size_type;
        typedef size_t                              difference_type;
        
        p_diagonal_matrix(size_t rows = 0, const T& init = T());
        size_t num_rows() const;
        size_t num_cols() const;
        p_diagonal_matrix<T>& operator = (const p_diagonal_matrix<T>& rhs);
        const T& operator[](size_t i) const;
        T & operator[](size_t i); 
        const T& operator()(size_t i, size_t j) const;
        T & operator()(size_t i, size_t j);
        void remove_rows(size_t i, size_t k = 1);
        void remove_cols(size_t j, size_t k = 1);
        void resize(size_t rows, size_t cols, T v = T());
	template< class T1 > 
        friend std::ostream & operator <<(std::ostream& os, const p_diagonal_matrix<T1>& m);
        const container& get_data() const; 
        container & get_data();    
        std::size_t size() const;
        void sqrt();
   private:
        container data_;
    };
}

#include "p_dense_matrix/p_diagonal_matrix.hpp"
#endif
