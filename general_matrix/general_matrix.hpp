#ifndef __ALPS_GENERAL_MATRIX_HPP__
#define __ALPS_GENERAL_MATRIX_HPP__

#include "detail/blasmacros.h"
#include "strided_iterator.hpp"
#include "vector.hpp"
#include "detail/general_matrix_adaptor.hpp"

#include <ostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>

//
// general_matrix template class
//
namespace blas {

    template <typename T, typename MemoryBlock = std::vector<T> >
    class general_matrix {
    public:
        // typedefs required for a std::container concept
        typedef T                       value_type;
        typedef T&                      reference;
        typedef T const&                const_reference;
        typedef std::size_t             size_type;
        typedef std::ptrdiff_t          difference_type;

        // TODO
        // for compliance with an std::container one would also need
        // -operators == != < > <= >=
        // -size()
        // -typedefs iterator, const_iterator
        // is it useful to implement this?

        // typedefs for matrix specific iterators
        // row_iterator: iterates over the rows of a specific column
        typedef strided_iterator<general_matrix,value_type>                  row_element_iterator;
        typedef strided_iterator<const general_matrix,const value_type>      const_row_element_iterator;
        // column_iterator: iterates over the columns of a specific row
        typedef value_type*                                                  column_element_iterator;
        typedef value_type const*                                            const_column_element_iterator;



        general_matrix(size_type size1 = 0, size_type size2 = 0, T init_value = T(0) )
        : size1_(size1), size2_(size2), reserved_size1_(size1), values_(size1*size2, init_value)
        {
        }

        general_matrix(general_matrix const& m)
        : size1_(m.size1_), size2_(m.size2_), reserved_size1_(m.size1_), values_()
        {
            // If the size of the matrix corresponds to the allocated size of the matrix...
            if(!m.is_shrinkable())
            {
                values_ = m.values_;
            }
            else
            {
                // copy only a shrinked to size version of the original matrix
                values_.reserve(m.size1_*m.size2_);
                for(size_type j=0; j < m.size2_; ++j)
                {
                    std::pair<const_column_element_iterator,const_column_element_iterator> range(m.column(j));
                    values_.insert(values_.end(), range.first, range.second);
                }
            }
        }

        template <typename OtherMemoryBlock>
        general_matrix(general_matrix<T,OtherMemoryBlock> const& m)
        : size1_(m.size1_), size2_(m.size2_), reserved_size1_(m.size1_), values_()
        {
            // If the size of the matrix corresponds to the allocated size of the matrix...
            if(!m.is_shrinkable())
            {
                values_.insert(values_.end(), m.values_.begin(), m.values_.end());
            }
            else
            {
                // copy only a shrinked to size version of the original matrix
                values_.reserve(m.size1_*m.size2_);
                for(size_type j=0; j < m.size2_; ++j)
                {
                    std::pair<typename general_matrix<T,OtherMemoryBlock>::const_column_element_iterator,
                              typename general_matrix<T,OtherMemoryBlock>::const_column_element_iterator
                             > range(m.column(j));
                    values_.insert(values_.end(), range.first, range.second);
                }
            }
        }
        
        void swap(general_matrix & r)
        {
            std::swap(values_, r.values_);
            std::swap(size1_, r.size1_);
            std::swap(size2_, r.size2_);
            std::swap(reserved_size1_,r.reserved_size1_);
        }
        
        friend void swap(general_matrix & x, general_matrix & y)
        {
            x.swap(y);
        }
       
        general_matrix& operator = (general_matrix rhs)
        {
            // swap(rhs, *this); // anyone have any idea why this doesn't work?
            this->swap(rhs);
            return *this;
        }

        inline value_type &operator()(const size_type i, const size_type j)
        {
            assert((i < size1_) && (j < size2_));
            return values_[i+j*reserved_size1_];
        }
        
        inline const value_type &operator()(const size_type i, const size_type j) const 
        {
            assert((i < size1_) && (j < size2_));
            return values_[i+j*reserved_size1_];
        }

        inline const bool empty() const
        {
            return (size1_ == 0 || size2_ == 0);
        }

        inline const size_type num_rows() const
        {
            return size1_;
        }

        inline const size_type num_columns() const
        {
            return size2_;
        }
        
        //TODO: shall these two functions be kept for compatibility or can we drop them?
        inline const size_type size1() const 
        {
            return size1_;
        }
  
        inline const size_type size2() const
        { 
            return size2_;
        }       
        
        inline const difference_type stride1() const
        {
            return 1;
        }
        
        inline const difference_type stride2() const
        {
            return reserved_size1_;
        }

        inline void resize(size_type size1, size_type size2, T init_value = T())
        {
           // Resizes the matrix to the size1 and size2 and allocates enlarges the vector if needed
           // If the new size for any dimension is smaller only elements outside the new size will be deleted.
           // If the new size is larger for any dimension the new elements will be initialized by the init_value.
           // All other elements will keep their value.

            // TODO: Over-resize matrix to 1.4 or 2 times the requested size
            if(size1 <= reserved_size1_)
            {
                //TODO What happens if resize or fill throw an exception?
                values_.resize(reserved_size1_*size2,init_value);

                if(size1 > size1_)
                {
                    // Reset all new elements which are in already reserved rows of already existing columns to init_value
                    // For all elements of new columns this is already done by values_.resize()
                    for(size_type j=0; j < size2; ++j)
                    {
                        std::fill(values_.begin()+j*reserved_size1_ + size1_, values_.begin()+j*reserved_size1_ + size1, init_value);
                    }
                }

            }
            else // size1 > reserved_size1_
            {
                // This is exception safe: If an exception is thrown, values_ and tmp won't get swapped.

                MemoryBlock tmp(size1*size2,init_value);
                for(size_type j=0; j < size2_; ++j)
                {
                    // Copy column by column
                    std::copy( values_.begin()+j*reserved_size1_, values_.begin()+j*reserved_size1_+size1_, tmp.begin()+j*size1);
                }
                std::swap(values_,tmp);
                reserved_size1_ = size1;
            }
            size1_=size1;
            size2_=size2;
        }
        
        inline void reserve(size_type size1, size_type size2)
        {
            if(size1*size2 > values_.capacity() )
            {
                    values_.reserve(size1*size2);
            }
            if(size1 > reserved_size1_)
            {
                MemoryBlock tmp(size1*size2);
                for(size_type j=0; j < size2_; ++j)
                {
                    // Copy column by column
                    std::copy( values_.begin()+j*reserved_size1_, values_.begin()+j*reserved_size1_+size1_, tmp.begin()+j*size1);
                }
                std::swap(values_,tmp);
                reserved_size1_ = size1;
            }
        }

        std::pair<size_type,size_type> capacity() const
        {
            assert( values_.capacity() % reserved_size1_ == 0 );
            // Evaluate the maximal number of columns (with size reserved_size1_) that the underlying vector could hold.
            // If the constructor, resize() and reserve() of std::vector would guarantee to allocate 
            // the requested amount of memory exactly
            // values_.capacity() % reserved_size1_ == 0 should hold.
            // However these functions guarantee to allocate _at least_ the requested amount.
            size_type reserved_size2_ = values_.capacity() - (values_.capacity() % reserved_size1_) / reserved_size1_;
            return std::pair<size_type,size_type>( reserved_size1_, reserved_size2_ );
        }
        
        bool is_shrinkable() const
        {
            // This assertion should actually never fail
            assert( reserved_size1_*size2_ == values_.size() );
            if(size1_ < reserved_size1_) return true;
            else return false;
        }

        void clear()
        {
            // Clear the values vector and ensure the reserved size stays the way it was
            values_.clear();
            values_.resize(reserved_size1_*size2_);
            size1_ = 0;
            size2_ = 0;
        }

        std::pair<row_element_iterator,row_element_iterator> row(size_type row = 0)
        {
            return std::make_pair( row_element_iterator(&values_[row],reserved_size1_), row_element_iterator(&values_[row+reserved_size1_*size2_], reserved_size1_) );
        }

        std::pair<const_row_element_iterator,const_row_element_iterator> row(size_type row = 0) const
        {
            return std::make_pair( const_row_element_iterator(&values_[row],reserved_size1_), const_row_element_iterator(&values_[row+reserved_size1_*size2_], reserved_size1_) );
        }

        std::pair<column_element_iterator,column_element_iterator> column(size_type col = 0 )
        {
            return std::make_pair( column_element_iterator(&values_[col*reserved_size1_]), column_element_iterator(&values_[col*reserved_size1_+size1_]) );
        }
        
        std::pair<const_column_element_iterator,const_column_element_iterator> column(size_type col = 0 ) const
        {
            return std::make_pair( const_column_element_iterator(&values_[col*reserved_size1_]), const_column_element_iterator(&values_[col*reserved_size1_+size1_]) );
        }

        template <typename InputIterator>
        void append_column(std::pair<InputIterator,InputIterator> const& range)
        {
            assert( std::distance(range.first, range.second) == size1_ );
            size_type insert_position = size2_;
            resize(size1_,size2_+1);    // This call modifies size2_ !
            std::copy( range.first, range.second, column(insert_position).first );
        }

        template <typename InputIterator>
        void append_row(std::pair<InputIterator,InputIterator> const& range)
        {
            assert( std::distance(range.first, range.second) == size2_ );
            size_type insert_position = size1_;
            resize(size1_+1,size2_);   // This call modifies size1_ !
            std::copy( range.first, range.second, row(insert_position).first );
        }

        template <typename InputIterator>
        void insert_row(size_type i, std::pair<InputIterator,InputIterator> const& range)
        {
            assert( i <= size1_ );
            assert( std::distance(range.first, range.second) == size2_ );

            // Append the row
            append_row(range);

            // Move the row through the matrix to the right possition
            for(size_type k=size1_-1; k>i; ++k)
            {
                swap_rows(k,k-1);
            }
        }

        template <typename InputIterator>
        void insert_column(size_type j, std::pair<InputIterator,InputIterator> const& range)
        {
            assert( j <= size2_);
            assert( std::distance(range.first, range.second) == size1_ );
            
            // Append the column
            append_column(range);

            // Move the column through the matrix to the right possition
            for(size_type k=size2_-1; k>j; ++k)
            {
                swap_columns(k,k-1);
            }

        }
        
        void swap_rows(size_type i1, size_type i2)
        {
            assert( i1 < size1_ && i2 < size1_ );
            std::pair<row_element_iterator, row_element_iterator> range( row(i1) );
            std::swap_ranges( range.first, range.second, row(i2).second );
        }

        void swap_columns(size_type j1, size_type j2)
        {
            assert( j1 < size2_ && j2 < size2_ );
            std::pair<column_element_iterator, column_element_iterator> range( column(j1) );
            std::swap_ranges(range.first, range.second, column(j2).first );
        }

        bool operator == (general_matrix const& rhs) const
        {
            if(size1_ != rhs.size1_ || size2_ != rhs.size2_) return false;
            // TODO: reimplement - this is just a quick ugly implementation
            for(size_type j=0; j < size2_; ++j)
                for(size_type i=0; i< size1_; ++i)
                    if(operator()(i,j) != rhs(i,j)) return false;
            return true;
        }

        general_matrix<T,MemoryBlock>& operator += (general_matrix const& rhs) 
        {
            using blas::plus_assign;
            plus_assign(*this,rhs);
            return *this;
        }
        
        general_matrix<T,MemoryBlock>& operator -= (general_matrix const& rhs) 
        {
            using blas::minus_assign;
            minus_assign(*this,rhs);
            return *this;
        }
        
        general_matrix<T,MemoryBlock>& operator *= (T const& t)
        {
            using blas::multiplies_assign;
            multiplies_assign(*this, t);
            return *this;
        }

        // Default implementations
        void plus_assign(general_matrix const& rhs)
        {
            assert((rhs.size1() == size1_) && (rhs.size2() == size2_));
            if(!(this->is_shrinkable() || rhs.is_shrinkable()) )
            {
                std::transform(this->values_.begin(),this->values_.end(),rhs.values_.begin(),this->values_.begin(), std::plus<T>());
            }
            else
            {
                // Do the operation column by column
                for(size_type j=0; j < size2_; ++j)
                {
                    std::pair<column_element_iterator,column_element_iterator> range(column(j));
                    std::transform( range.first, range.second, rhs.column(j).first, range.first, std::plus<T>());
                }
            }
        }

        void minus_assign(general_matrix const& rhs)
        {
            assert((rhs.size1() == size1_) && (rhs.size2() == size2_));
            if(!(this->is_shrinkable() || rhs.is_shrinkable()) )
            {
                std::transform(this->values_.begin(),this->values_.end(),rhs.values_.begin(),this->values_.begin(), std::minus<T>());
            }
            else
            {
                // Do the operation column by column
                for(size_type j=0; j < size2_; ++j)
                {
                    std::pair<column_element_iterator,column_element_iterator> range(column(j));
                    std::transform( range.first, range.second, rhs.column(j).first, range.first, std::minus<T>());
                }
            }
        }
        
        void multiplies_assign (T const& t)
        {
            if(!(is_shrinkable()) )
            {
                std::transform(values_.begin(),values_.end(),values_.begin(), bind1st(std::multiplies<T>(),t));
            }
            else
            {
                // Do the operation column by column
                for(size_type j=0; j < size2_; ++j)
                {
                    std::pair<column_element_iterator,column_element_iterator> range(column(j));
                    std::transform(range.first, range.second, range.first, bind1st(std::multiplies<T>(),t));
                }
            }
        }
        
        general_matrix<T,MemoryBlock>& operator *= (general_matrix const& rhs) 
        {
            // It's not common to implement a *= operator in terms of a * operator,
            // but a temporary object has to be created to store the result anyway
            // so it seems reasonable.
            general_matrix tmp = (*this) * rhs;
            std::swap(tmp,*this);
            return *this;
        }

    private:
        template <typename OtherT,typename OtherMemoryBlock>
        friend class general_matrix;

//        friend class boost::numeric::bindings::detail::adaptor<general_matrix<T,MemoryBlock>,const general_matrix<T,MemoryBlock>, void>;
//        friend class boost::numeric::bindings::detail::adaptor<general_matrix<T,MemoryBlock>,general_matrix<T,MemoryBlock>, void>;


        size_type size1_;
        size_type size2_;
        size_type reserved_size1_;
        // "reserved_size2_" is done automatically by underlying std::vector (see vector.reserve(), vector.capacity() )
        
        MemoryBlock values_;
    };
} // namespace blas


//
// Function hooks
//

namespace blas {

    template <typename T, typename MemoryBlock>
    const general_matrix<T,MemoryBlock> matrix_matrix_multiply(general_matrix<T,MemoryBlock> const& lhs, general_matrix<T,MemoryBlock> const& rhs)
    {
        assert( lhs.num_columns() == rhs.num_rows() );

        // Simple matrix matrix multiplication
        general_matrix<T,MemoryBlock> result(lhs.num_rows(),rhs.num_columns());
        for(std::size_t i=0; i < lhs.num_rows(); ++i)
        {
            for(std::size_t j=0; j<rhs.num_columns(); ++j)
            {
                for(std::size_t k=0; k<lhs.num_columns(); ++k)
                {
                        result(i,j) += lhs(i,k) * rhs(k,j);
                }
            }
        }
        return result;
    } 
    
    template <typename T,typename MemoryBlock>
    void plus_assign(general_matrix<T,MemoryBlock>& m, general_matrix<T,MemoryBlock> const& rhs)
    {
        m.plus_assign(rhs);
    }

    template <typename T, typename MemoryBlock>
    void minus_assign(general_matrix<T,MemoryBlock>& m, general_matrix<T,MemoryBlock> const& rhs)
    {
        m.minus_assign(rhs);
    }

    template <typename T, typename MemoryBlock>
    void multiplies_assign(general_matrix<T,MemoryBlock>& m, T const& t)
    {
        m.multiplies_assign(t);
    }
}

//
// Free general matrix functions
//
namespace blas {
/*
    template<typename MatrixType, class DoubleVector>
    void svd(MatrixType M, MatrixType& U, MatrixType& V, DoubleVector & S)
    {
        std::size_t K = std::min(M.num_rows(), M.num_columns());
        U.resize(M.num_rows(), K);
        V.resize(K, M.num_columns());
        S.resize(K);
        
        boost::numeric::bindings::lapack::gesdd('S', M, S, U, V);
    }
*/

    // 
    template <typename MatrixType>
    MatrixType transpose(MatrixType const& m) 
    {
        // TODO: perhaps this could return a proxy object
        MatrixType tmp(m.size2(),m.size1());
        for(std::size_t i=0;i<m.size1();++i){
            for(std::size_t j=0;j<m.size2();++j){
                tmp(j,i) = m()(i,j);
            }
        }
        return tmp;
    }
    
    template <typename MatrixType>
    const typename MatrixType::value_type trace(MatrixType const& m)
    {
        assert(m.size1() == m.size2());
        typename MatrixType::value_type tr(0);
        for(std::size_t i=0; i<m.size1(); ++i) tr+=m(i,i);
        return tr;
    }
    
    template <typename T, typename MemoryBlock>
    const general_matrix<T,MemoryBlock> operator + (general_matrix<T,MemoryBlock> a, general_matrix<T,MemoryBlock> const& b)
    {
        a += b;
        return a;
    }
    
    template <typename T, typename MemoryBlock>
    const general_matrix<T,MemoryBlock> operator - (general_matrix<T,MemoryBlock> a, general_matrix<T,MemoryBlock> const& b)
    {
        a -= b;
        return a;
    }

    template<typename T, typename MemoryBlock>
    const vector<T,MemoryBlock> operator * (general_matrix<T,MemoryBlock> const& m, vector<T,MemoryBlock> const& v)
    {
        assert( m.size2() == v.size() );
        vector<T,MemoryBlock> result(m.size1());
        // Simple Matrix * Vector
        for(typename general_matrix<T,MemoryBlock>::size_type i = 0; i < m.size1(); ++i)
        {
            for(typename general_matrix<T,MemoryBlock>::size_type j=0; j <m.size2(); ++j)
            {
                result(i) = m(i,j) * v(j);
            }
        }
        return result;
    }
   
    // TODO: adj(Vector) * Matrix, where adj is a proxy object


    template<typename T,typename MemoryBlock>
    const general_matrix<T,MemoryBlock> operator * (general_matrix<T,MemoryBlock> m, T const& t)
    {
        return m*=t;
    }
    
    template<typename T,typename MemoryBlock>
    const general_matrix<T,MemoryBlock> operator * (T const& t, general_matrix<T,MemoryBlock> m)
    {
        return m*=t;
    }

    template<typename T, typename MemoryBlock>
    const general_matrix<T,MemoryBlock> operator * (general_matrix<T,MemoryBlock> const& m1, general_matrix<T,MemoryBlock> const& m2)
    {
        return matrix_matrix_multiply(m1,m2);
    }

    template<typename T,typename MemoryBlock>
    void gemm(general_matrix<T,MemoryBlock> const & A, general_matrix<T,MemoryBlock> const & B, general_matrix<T,MemoryBlock> & C)
    {
        C = matrix_matrix_multiply(A, B);
    }

    template <typename T, typename MemoryBlock>
    std::ostream& operator << (std::ostream& o, general_matrix<T,MemoryBlock> const& m)
    {
        for(typename general_matrix<T,MemoryBlock>::size_type i=0; i< m.num_rows(); ++i)
        {
            for(typename general_matrix<T,MemoryBlock>::size_type j=0; j < m.num_columns(); ++j)
                o<<m(i,j)<<" ";
            o<<std::endl;
        }
        return o;
    }
} // namespace blas

#endif //__ALPS_GENERAL_MATRIX_HPP__
