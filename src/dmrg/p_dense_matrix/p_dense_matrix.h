#ifndef __ALPS_DENSE_MATRIX_HPP__
#define __ALPS_DENSE_MATRIX_HPP__

#include "ambient/ambient.hpp"

#include "p_dense_matrix/iterators.hpp"
#include "p_dense_matrix/vector.hpp"

#include "utils/zout.hpp"
#include "utils/function_objects.h"

#include "p_dense_matrix/p_diagonal_matrix.h"

#include <boost/lambda/lambda.hpp>
#include <boost/typeof/typeof.hpp>
#include <ostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
//#include <alps/ngs/hdf5/deprecated.hpp>
#endif

namespace blas {

    template <typename T, ambient::policy P>
    class p_dense_matrix : public ambient::livelong<p_dense_matrix<T,P>, P>
    {
    public:
        typedef p_dense_matrix<T, ambient::REPLICA> replica;
        typedef ambient::livelong<p_dense_matrix<T,P>, P> livelong;
        friend class p_dense_matrix<T, ambient::ANY>;
        friend class p_dense_matrix<T, ambient::REPLICA>;
        ~p_dense_matrix();
// AMBIENT PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        p_dense_matrix(ambient::void_pt* p);              // proxy object construction
// AMBIENT PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        typedef T         value_type;                     // The type T of the elements of the matrix
        typedef size_t    size_type;                      // Unsigned integer type that represents the dimensions of the matrix
        typedef ptrdiff_t difference_type;                // Signed integer type to represent the distance of two elements in the memory

        // typedefs for matrix specific iterators
        typedef strided_iterator<p_dense_matrix,value_type>
            row_element_iterator;                         // Iterator to iterate through the elements of a row of the matrix
        typedef strided_iterator<const p_dense_matrix,const value_type>
            const_row_element_iterator;                   // Const version of row_element_iterator
        typedef value_type*
            column_element_iterator;                      // Iterator to iterate through the elements of a cols of the matrix
        typedef value_type const*
            const_column_element_iterator;                // Const version of column_element_iterator       
        typedef matrix_element_iterator<p_dense_matrix,value_type>
            element_iterator;                             // Iterator to iterate through all elements of the matrix (REALLY SLOW! USE row_-/column_iterators INSTEAD!)
        typedef matrix_element_iterator<const p_dense_matrix,const value_type>
            const_element_iterator;                       // Const version of element_iterator (REALLY SLOW! USE row_-/column_iterators INSTEAD!)

        p_dense_matrix(size_type rows, size_type cols, T init_value);
        template<typename TM, ambient::policy PM>
        p_dense_matrix(p_dense_matrix<TM,PM> const& m);
  
        void swap(p_dense_matrix & r);
        friend void swap(p_dense_matrix & x, p_dense_matrix & y){ x.swap(y); }

        inline bool empty() const;
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        inline difference_type stride1() const;
        inline difference_type stride2() const;
        inline difference_type get_lda() const;
        inline difference_type get_sda() const;
        void resize(size_type rows, size_type cols);
        void remove_rows(size_type i, size_type k);
        void remove_cols(size_type j, size_type k);
        void clear();
       
        void inplace_conjugate();

        inline value_type& get(size_type i, size_type j) const;

        template<ambient::policy PR> operator p_dense_matrix<T,PR> ();
        template<ambient::policy PR> p_dense_matrix<T,P>& operator = (p_dense_matrix<T,PR>& rhs);
        p_dense_matrix<T,P>& operator = (const p_dense_matrix<T>& rhs);
        inline value_type& operator()(size_type i, size_type j);
        inline const value_type& operator()(size_type i, size_type j) const;
        p_dense_matrix<T,P>& operator += (const p_dense_matrix& rhs); 
        p_dense_matrix<T,P>& operator -= (const p_dense_matrix& rhs);
        template <typename T2>
        p_dense_matrix<T,P>& operator *= (const T2& t);
        template <typename T2>
        p_dense_matrix<T,P>& operator /= (const T2& t);

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
            return std::make_pair( element_iterator(this->self,0,0), element_iterator(this->self,0,num_cols()) );
        }
        std::pair<element_iterator,element_iterator> elements() const
        {
            return std::make_pair( const_element_iterator(this->self,0,0), const_element_iterator(this->self,0,num_cols() ) );
        }

        element_iterator begin()
        {
            return element_iterator(this->self,0,0);
        }

        const_element_iterator begin() const 
        {
            return const_element_iterator(this->self,0,0);
        }

        element_iterator end()
        {
            return element_iterator(this->self,num_rows(),num_cols());
        } 

        const_element_iterator end() const
        {
            return const_element_iterator(this->self,num_rows(),num_cols());
        } 
/*
        element_iterator push_back()
        {
          assert(false);
          return element_iterator;
        }

        const_element_iterator push_back() const
        {
          assert(false);
          return element_iterator;
        }
*/

#ifdef HAVE_ALPS_HDF5
        void serialize(alps::hdf5::iarchive & ar);
        void serialize(alps::hdf5::oarchive & ar) const;
#endif
    private:
        size_type rows;
        size_type cols;
        size_type lda; // leading dimension array
        size_type sda; // subleading dimension array
        T* data;
    };

    template<typename T, ambient::policy P>
    struct associated_diagonal_matrix< p_dense_matrix<T,P> >
    {
        typedef p_diagonal_matrix<T> type;
    };
    
} // namespace blas

namespace blas {

    template<class FullMatrixClass>
    struct associated_vector { };

    template<typename T, ambient::policy P>
    struct associated_vector<p_dense_matrix<T,P> >
    {
        typedef p_diagonal_matrix<T> type;
    };
} // namespace blas



#include "p_dense_matrix/p_dense_matrix.hpp"
#endif //__ALPS_DENSE_MATRIX_HPP__
