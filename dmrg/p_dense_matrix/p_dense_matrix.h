#ifndef __ALPS_DENSE_MATRIX_HPP__
#define __ALPS_DENSE_MATRIX_HPP__

#include "ambient/ambient.h"

#include "p_dense_matrix/iterators.hpp"
#include "p_dense_matrix/vector.hpp"

#include "utils/function_objects.h"

#include "dense_matrix/diagonal_matrix.h"

#include <boost/scoped_ptr.hpp>
#include <ostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

namespace blas {
    template <typename T>
    class p_dense_matrix {
    public:
        typedef T                       value_type;       // The type T of the elements of the matrix
        typedef std::size_t             size_type;        // Unsigned integer type that represents the dimensions of the matrix
        typedef std::ptrdiff_t          difference_type;  // Signed integer type to represent the distance of two elements in the memory

        // typedefs for matrix specific iterators
        typedef strided_iterator<p_dense_matrix,value_type>
            row_element_iterator;                         // Iterator to iterate through the elements of a row of the matrix
        typedef strided_iterator<const p_dense_matrix,const value_type>
            const_row_element_iterator;                   // Const version of row_element_iterator
        typedef value_type*
            column_element_iterator;                      // Iterator to iterate through the elements of a columns of the matrix
        typedef value_type const*
            const_column_element_iterator;                // Const version of column_element_iterator       
        typedef matrix_element_iterator<p_dense_matrix,value_type>
            element_iterator;                             // Iterator to iterate through all elements of the matrix (REALLY SLOW! USE row_-/column_iterators INSTEAD!)
        typedef matrix_element_iterator<const p_dense_matrix,const value_type>
            const_element_iterator;                       // Const version of element_iterator (REALLY SLOW! USE row_-/column_iterators INSTEAD!)

        p_dense_matrix(size_type rows, size_type columns, T init_value);
        p_dense_matrix(p_dense_matrix const& m);
        void swap(p_dense_matrix & r);
        friend void swap(p_dense_matrix & x, p_dense_matrix & y)
        {
            x.swap(y);
        }

        p_dense_matrix& operator = (p_dense_matrix rhs);
        inline value_type& operator()(const size_type i, const size_type j);
        inline value_type const& operator()(const size_type i, const size_type j) const;
        p_dense_matrix<T>& operator += (p_dense_matrix const& rhs); 
        p_dense_matrix<T>& operator -= (p_dense_matrix const& rhs);
        
        template <typename T2>
        p_dense_matrix<T>& operator *= (T2 const& t);
        
        template <typename T2>
        p_dense_matrix<T>& operator /= (T2 const& t);

        inline const bool empty() const;
        inline const size_type num_rows() const;
        inline const size_type num_columns() const;
        inline const difference_type stride1() const;
        inline const difference_type stride2() const;
        void reserve(size_type rows, size_type cols);
        void resize(size_type rows, size_type cols);
        void clear();

        std::pair<row_element_iterator,row_element_iterator> row(size_type row = 0)
        {
            return std::make_pair( row_element_iterator(&data[row],lda), row_element_iterator(&data[row+lda*cols], lda) );
        }
        std::pair<const_row_element_iterator,const_row_element_iterator> row(size_type row = 0) const
        {
            return std::make_pair( const_row_element_iterator(&data[row],lda), const_row_element_iterator(&data[row+lda*cols], lda) );
        }
        std::pair<column_element_iterator,column_element_iterator> column(size_type col = 0 )
        {
            return std::make_pair( column_element_iterator(&data[col*lda]), column_element_iterator(&data[col*lda+rows]) );
        }
        std::pair<const_column_element_iterator,const_column_element_iterator> column(size_type col = 0 ) const
        {
            return std::make_pair( const_column_element_iterator(&data[col*lda]), const_column_element_iterator(&data[col*lda+rows]) );
        }
        std::pair<element_iterator,element_iterator> elements()
        {
            return std::make_pair( element_iterator(this,0,0), element_iterator(this,0, num_columns()) );
        }
        std::pair<element_iterator,element_iterator> elements() const
        {
            return std::make_pair( const_element_iterator(this,0,0), const_element_iterator(this,0,num_columns() ) );
        }

        void remove_rows(size_type i, difference_type k);
        void remove_columns(size_type j, difference_type k);
        void plus_assign(p_dense_matrix const& rhs);
        void minus_assign(p_dense_matrix const& rhs);

        template <typename T2>
        void multiplies_assign (T2 const& t);
       
        void inplace_conjugate();
#ifdef HAVE_ALPS_HDF5
        void serialize(alps::hdf5::iarchive & ar);
        void serialize(alps::hdf5::oarchive & ar) const;
#endif
    private:
        size_type rows;
        size_type cols;
        size_type lda; // leading dimension array
        size_type sda; // subleading dimension array
        boost::scoped_ptr<T> data_scope;
        T* data;
    };
    
    template<typename T>
    struct associated_diagonal_matrix< p_dense_matrix<T> >
    {
        typedef diagonal_matrix<T> type;
    };
    
} // namespace blas


namespace blas {
    
    template <typename T>
    const p_dense_matrix<T> operator + (p_dense_matrix<T> a, p_dense_matrix<T> const& b);
    
    template <typename T>
    const p_dense_matrix<T> operator - (p_dense_matrix<T> a, p_dense_matrix<T> const& b);

    template<typename T, typename T2>
    const p_dense_matrix<T> operator * (p_dense_matrix<T> m, T2 const& t);
    
    template<typename T, typename T2>
    const p_dense_matrix<T> operator * (T2 const& t, p_dense_matrix<T> m);

    template<typename T>
    const p_dense_matrix<T> operator * (p_dense_matrix<T> const& m1, p_dense_matrix<T> const& m2);

    template <typename T>
    std::ostream& operator << (std::ostream& o, p_dense_matrix<T> const& m);

} // namespace blas

#include "p_dense_matrix/p_dense_matrix.hpp"
#endif //__ALPS_DENSE_MATRIX_HPP__
