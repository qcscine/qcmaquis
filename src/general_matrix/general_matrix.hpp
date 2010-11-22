#ifndef __ALPS_GENERAL_MATRIX_HPP__
#define __ALPS_GENERAL_MATRIX_HPP__

#include "strided_iterator.hpp"
#include "vector.hpp"
#include "detail/general_matrix_adaptor.hpp"

#include <boost/lambda/lambda.hpp>
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

        // for compliance with an std::container one would also need
        // -operators == != < > <= >=
        // -size()
        // -typedefs iterator, const_iterator


        // typedefs for matrix specific iterators
        // row_iterator: iterates over the rows of a specific column
        typedef strided_iterator<general_matrix,value_type>                  row_element_iterator;
        typedef strided_iterator<const general_matrix,const value_type>      const_row_element_iterator;
        // column_iterator: iterates over the columns of a specific row
        typedef value_type*                                                  column_element_iterator;
        typedef value_type const*                                            const_column_element_iterator;



        general_matrix(size_type size1 = 0, size_type size2 = 0, T init_value = T() )
        : size1_(size1), size2_(size2), reserved_size1_(size1), values_(size1*size2, init_value)
        {
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

        void resize(size_type size1, size_type size2, T const& init_value = T())
        {
           // Resizes the matrix to the size1 and size2 and enlarges the
           // MemoryBlock if needed. If the new size for any dimension is
           // smaller only elements outside the new size will be deleted.
           // If the new size is larger for any dimension the new elements
           // will be initialized by the init_value.
           // All other elements will keep their value.
            
            // Exception behaviour:
            // As long as the assignment and copy operation of the T values don't throw an exception,
            // any exception will leave the matrix unchanged.
            // (Assuming the same behaviour of the underlying MemoryBlock. This is true for std::vector.)

            // Do we need more space? Reserve more space if needed!
            //
            // If the memory is reallocated using reserve
            // we just have to fill the new columns with the init_value
            // (->after this if statement),
            // since reserve fills all elements in the range between size1_
            // and reserved_size1_ of each EXISTING column with init_value
            // by using values_.resize()
            if(!automatic_reserve(size1,size2,init_value))
            {
                if(size1 > size1_)
                {
                    // Reset all "new" elements which are in already reserved
                    // rows of already existing columns to init_value
                    // For all elements of new columns this is already done by
                    // values_.resize() (->after this if statement)
                    size_type num_of_cols = std::min(size2, size2_);
                    for(size_type j=0; j < num_of_cols; ++j)
                        std::fill(
                                values_.begin()+j*reserved_size1_ + size1_,
                                values_.begin()+j*reserved_size1_ + size1,
                                init_value
                                );
                }
            }
            values_.resize(reserved_size1_*size2, init_value);
            size1_=size1;
            size2_=size2;
        }
        
        void reserve(size_type size1, size_type size2, T const& init_value = T())
        {
            // The init_value may seem a little weird in a reserve method,
            // but one has to initialize all matrix elements in the
            // reserved_size1_ range of each column, due to the 1d-structure
            // of the underlying MemoryBlock (e.g. std::vector)

            // Ignore values that would shrink the matrix
            size2 = std::max(size2, size2_);
            size1 = std::max(size1, reserved_size1_);
           
            // Is change of structure or size of the MemoryBlock necessary?
            if(size1 > reserved_size1_ || size1*size2 > values_.capacity() )
            {
                MemoryBlock tmp;
                tmp.reserve(size1*size2);
                // Copy column by column
                for(size_type j=0; j < size2_; ++j)
                {
                    std::pair<column_element_iterator, column_element_iterator> range(column(j));
                    // Copy the elements from the current MemoryBlock
                    tmp.insert(tmp.end(),range.first,range.second);
                    // and fill the rest with the init_value
                    tmp.insert(tmp.end(),size1-size1_,init_value);
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
        void append_columns(std::pair<InputIterator,InputIterator> const& range, difference_type k=1)
        {
            assert( std::distance(range.first, range.second) == k*size1_ );
            // Reserve more space if needed
            automatic_reserve(size1_,size2_+k);
            // Append column by column
            for(difference_type l=0; l<k; ++l)
            {
                values_.insert(values_.end(), range.first+(l*size1_), range.first+((l+1)*size1_) );
                // Fill the space reserved for new rows
                values_.insert(values_.end(), reserved_size1_-size1_, T());
            }
            size2_ += k;
        }

        template <typename InputIterator>
        void append_rows(std::pair<InputIterator,InputIterator> const& range, difference_type k =1)
        {
            assert( std::distance(range.first, range.second) == k*size2_ );
            // Reserve more space if needed
            automatic_reserve(size1_+k, size2_);
            // The elements do already exists due to reserve, so we can just use (copy to) them.
            for(difference_type l=0; l<k; ++l)
                std::copy( range.first+(l*size2_), range.first+((l+1)*size2_), row(size1_+l).first );
            size1_ += k;
        }

        template <typename InputIterator>
        void insert_rows(size_type i, std::pair<InputIterator,InputIterator> const& range, difference_type k = 1)
        {
            assert( i <= size1_ );
            assert( std::distance(range.first, range.second) == k*size2_ );

            // Append the row
            automatic_reserve(size1_+k,size2_);

            for(size_type j=0; j<size2_; ++j)
                std::copy_backward(&values_[reserved_size1_*j+i],&values_[reserved_size1_*j+size1_],&values_[reserved_size1_*j+size1_+k]);
            for(difference_type l=0; l<k; ++l)
                std::copy(range.first+l*size2_,range.first+(l+1)*size2_,row(i+l).first);
            size1_+=k;
        }

        template <typename InputIterator>
        void insert_columns(size_type j, std::pair<InputIterator,InputIterator> const& range, difference_type k = 1)
        {
            assert( j <= size2_);
            assert( std::distance(range.first, range.second) == k*size1_ );
            
            // Append the column
            automatic_reserve(size1_,size2_+k);

            // Move the column through the matrix to the right possition
            for(size_type h=size2_; h>j; --h)
                std::copy(&values_[reserved_size1_*(h-1)],&values_[reserved_size1_*(h-1)]+size1_,&values_[reserved_size1_*(h+k-1)]);
            for(difference_type l=0; l<k; ++l)
                std::copy(range.first+l*size1_,range.first+(l+1)*size1_,&values_[reserved_size1_*(j+l)]);
            size2_+=k;
        }

        void remove_rows(size_type i, difference_type k = 1)
        {
            assert( i+k <= size1_ );
            // for each column, copy the rows > i+k   k rows  up
            for(size_type j = 0; j < size2_; ++j)
                std::copy(&values_[reserved_size1_*j + i + k], &values_[reserved_size1_*j + size1_], &values_[reserved_size1_*j + i] );
            size1_ -= k;
        }

        void remove_columns(size_type j, difference_type k = 1)
        {
            assert( j+k <= size2_ );
            values_.erase(values_.begin()+(reserved_size1_*j), values_.begin()+(reserved_size1_*(j+k)) );
            size2_ -= k;
        }

        void swap_rows(size_type i1, size_type i2)
        {
            assert( i1 < size1_ && i2 < size1_ );
            std::pair<row_element_iterator, row_element_iterator> range( row(i1) );
            std::swap_ranges( range.first, range.second, row(i2).first );
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
            // TODO: reimplement or remove - this is just a quick ugly implementation
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
        
        template <typename T2>
        general_matrix<T,MemoryBlock>& operator *= (T2 const& t)
        {
            using blas::multiplies_assign;
            multiplies_assign(*this, t);
            return *this;
        }

        // Default implementations
        void plus_assign(general_matrix const& rhs)
        {
            assert((rhs.size1_ == size1_) && (rhs.size2_ == size2_));
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
            assert((rhs.size1_ == size1_) && (rhs.size2_ == size2_));
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

        template <typename T2>
        void multiplies_assign (T2 const& t)
        {
            if(!(is_shrinkable()) )
            {
                std::for_each(values_.begin(), values_.end(), boost::lambda::_1 *= t);
            }
            else
            {
                // Do the operation column by column
                for(size_type j=0; j < size2_; ++j)
                {
                    std::pair<column_element_iterator,column_element_iterator> range(column(j));
                    std::for_each(range.first, range.second, boost::lambda::_1 *= t);
                }
            }
        }
        
    private:
        template <typename OtherT,typename OtherMemoryBlock>
        friend class general_matrix;

        inline bool automatic_reserve(size_type size1, size_type size2, T const& init_value = T())
        {
            // Do we need to reserve more space in any dimension?
            if(size1 > reserved_size1_ || reserved_size1_*size2 > values_.capacity())
            {
                reserve(size1*3/2,size2*3/2,init_value);
                return true;
            }
            else
            {
                return false;
            }
        }

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

    template <typename T, typename MemoryBlock, typename T2>
    void multiplies_assign(general_matrix<T,MemoryBlock>& m, T2 const& t)
    {
        m.multiplies_assign(t);
    }
}

//
// Free general matrix functions
//
namespace blas {
    
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

    // TODO Return type deduction!
    template<typename T, typename MemoryBlock>
    const vector<T,MemoryBlock> operator * (general_matrix<T,MemoryBlock> const& m, vector<T,MemoryBlock> const& v)
    {
        assert( m.num_columns() == v.size() );
        vector<T,MemoryBlock> result(m.num_rows());
        // Simple Matrix * Vector
        for(typename general_matrix<T,MemoryBlock>::size_type i = 0; i < m.num_rows(); ++i)
        {
            for(typename general_matrix<T,MemoryBlock>::size_type j=0; j <m.num_columns(); ++j)
            {
                result(i) = m(i,j) * v(j);
            }
        }
        return result;
    }
   
    // TODO: adj(Vector) * Matrix, where adj is a proxy object

    template<typename T,typename MemoryBlock, typename T2>
    const general_matrix<T,MemoryBlock> operator * (general_matrix<T,MemoryBlock> m, T2 const& t)
    {
        return m*=t;
    }
    
    template<typename T,typename MemoryBlock, typename T2>
    const general_matrix<T,MemoryBlock> operator * (T2 const& t, general_matrix<T,MemoryBlock> m)
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
