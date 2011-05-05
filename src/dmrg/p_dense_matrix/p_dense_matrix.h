#ifndef __ALPS_P_DENSE_MATRIX_HPP__
#define __ALPS_P_DENSE_MATRIX_HPP__

#include "ambient/ambient.hpp"


#include "utils/zout.hpp"
#include "utils/function_objects.h"

#include "p_dense_matrix/iterators.hpp"
#include "p_dense_matrix/p_diagonal_matrix.h"

#include <ostream>

namespace blas {

    template <typename T, ambient::policy P>
    class p_dense_matrix : public ambient::livelong<p_dense_matrix<T,P>, P>
    {
    public:
        typedef p_dense_matrix<T, ambient::REPLICA> replica;
        typedef ambient::livelong<p_dense_matrix<T,P>, P> livelong;
        friend class p_dense_matrix<T, ambient::ANY>;
        friend class p_dense_matrix<T, ambient::REPLICA>;
        typedef T         value_type;                     // The type T of the elements of the matrix
        typedef size_t    size_type;                      // Unsigned integer type that represents the dimensions of the matrix
        typedef ptrdiff_t difference_type;                // Signed integer type to represent the distance of two elements in the memory

        // typedefs for matrix specific iterators // don't use iterators!! they are really slow!! reprogram with kernels instead.
        typedef p_matrix_element_iterator<p_dense_matrix,value_type>             row_element_iterator;          // Iterator to iterate through the elements of a row of the matrix
        typedef p_matrix_element_iterator<const p_dense_matrix,const value_type> const_row_element_iterator;    // Const version of row_element_iterator
        typedef p_matrix_element_iterator<p_dense_matrix,value_type>             column_element_iterator;       // Iterator to iterate through the elements of a cols of the matrix
        typedef p_matrix_element_iterator<const p_dense_matrix,const value_type> const_column_element_iterator; // Const version of column_element_iterator       
        typedef p_matrix_element_iterator<p_dense_matrix,value_type>             element_iterator;              // Iterator to iterate through all elements of the matrix
        typedef p_matrix_element_iterator<const p_dense_matrix,const value_type> const_element_iterator;        // Const version of element_iterator

       ~p_dense_matrix();
        p_dense_matrix(ambient::void_pt* p);              // proxy object construction
        p_dense_matrix(size_type rows, size_type cols, T init_value);
        template<typename TM, ambient::policy PM>
        p_dense_matrix(p_dense_matrix<TM,PM> const& m);
  
        void swap(p_dense_matrix & r);
        friend void swap(p_dense_matrix & x, p_dense_matrix & y){ x.swap(y); }

        inline bool empty() const;
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        void resize(size_type rows, size_type cols);
        void remove_rows(size_type i, size_type k);
        void remove_cols(size_type j, size_type k);
        void nullcut();
        void clear();
       
        void inplace_conjugate();

        inline value_type& get(size_type i, size_type j) const;

        template<ambient::policy PR>                      operator p_dense_matrix<T,PR> ();
        template<ambient::policy PR> p_dense_matrix<T,P>& operator = (p_dense_matrix<T,PR>& rhs);
        p_dense_matrix<T,P>&                              operator = (const p_dense_matrix<T>& rhs);
        inline value_type&                                operator ()(size_type i, size_type j);
        inline const value_type&                          operator()(size_type i, size_type j) const;
        p_dense_matrix<T,P>&                              operator += (const p_dense_matrix& rhs); 
        p_dense_matrix<T,P>&                              operator -= (const p_dense_matrix& rhs);
        template <typename T2> p_dense_matrix<T,P>&       operator *= (const T2& t);
        template <typename T2> p_dense_matrix<T,P>&       operator /= (const T2& t);

        std::pair<row_element_iterator,row_element_iterator> row(size_type row = 0){
            return std::make_pair( row_element_iterator(this->self,row,0), row_element_iterator(this->self,row,cols) );
        }
        std::pair<const_row_element_iterator,const_row_element_iterator> row(size_type row = 0) const {
            return std::make_pair( const_row_element_iterator(this->self,row,0), const_row_element_iterator(this->self,row,cols) );
        }
        std::pair<column_element_iterator,column_element_iterator> column(size_type col = 0 ){
            return std::make_pair( column_element_iterator(this->self, 0, col), column_element_iterator(this->self, rows, col) );
        }
        std::pair<const_column_element_iterator,const_column_element_iterator> column(size_type col = 0 ) const {
            return std::make_pair( const_column_element_iterator(this->self, 0, col), const_column_element_iterator(this->self, rows, col) );
        }
        std::pair<element_iterator,element_iterator> elements(){
            return std::make_pair( element_iterator(this->self,0,0), element_iterator(this->self,0,num_cols()) );
        }
        std::pair<element_iterator,element_iterator> elements() const {
            return std::make_pair( const_element_iterator(this->self,0,0), const_element_iterator(this->self,0,num_cols() ) );
        }

    private:
        size_type rows;
        size_type cols;
    };

    template<typename T, ambient::policy P>
    struct associated_diagonal_matrix< p_dense_matrix<T,P> >
    {
        typedef p_diagonal_matrix<T> type;
    };

    template<typename T, ambient::policy P>
    struct associated_vector<p_dense_matrix<T,P> >
    {
        typedef p_diagonal_matrix<T> type;
    };
} // namespace blas

#include "p_dense_matrix/p_dense_matrix.hpp"
#endif //__ALPS_DENSE_MATRIX_HPP__
