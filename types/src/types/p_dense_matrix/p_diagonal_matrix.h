#ifndef __MAQUIS_TYPES_P_DIAGONAL_MATRIX_H
#define __MAQUIS_TYPES_P_DIAGONAL_MATRIX_H

namespace maquis { namespace types {
   
    template<typename T>
    class p_dense_matrix;

    template<typename T>
    class p_diagonal_matrix
    {
    public:
        typedef typename maquis::types::p_dense_matrix<T> container;
        typedef typename container::difference_type difference_type;
        typedef typename container::value_type value_type;
        typedef typename container::size_type size_type;
        
        inline p_diagonal_matrix();
        inline p_diagonal_matrix(size_t rows, const value_type& init = value_type());
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline p_diagonal_matrix<T>& operator = (const p_diagonal_matrix<T>& rhs);
        inline const value_type& operator[](size_t i) const;
        inline value_type& operator[](size_t i); 
        inline const value_type& operator()(size_t i, size_t j) const;
        inline value_type& operator()(size_t i, size_t j);
        inline void remove_rows(size_t i, size_t k = 1);
        inline void remove_cols(size_t j, size_t k = 1);
        inline void resize(size_t rows, size_t cols);
        template< class T1 > 
        friend std::ostream & operator <<(std::ostream& os, const p_diagonal_matrix<T1>& m);
        inline const container& get_data() const; 
        inline container & get_data();    
        inline size_type size() const;
        inline void sqrt();
        inline void exp(const T& alfa = 1.);
   private:
        container data_;
    };

} } // namespace maquis::types

#include "types/p_dense_matrix/p_diagonal_matrix.hpp"
#endif
