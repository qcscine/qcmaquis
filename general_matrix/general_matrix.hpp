#ifndef __ALPS_GENERAL_MATRIX_HPP__
#define __ALPS_GENERAL_MATRIX_HPP__

#include "strided_iterator.hpp"
#include "matrix_element_iterator.hpp"
#include "vector.hpp"
#include "detail/general_matrix_adaptor.hpp"

#include "../util/function_objects.h"

#include "diagonal_matrix.h"

#include <boost/lambda/lambda.hpp>
#include <boost/typeof/typeof.hpp>
#include <ostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>

namespace blas {
    /** A matrix template class
      *
      * The general_matrix class is a matrix which can take any type T
      * @param T the type of the elements to be stored in the matrix
      * @param MemoryBlock the underlying (continous) Memory structure
      */
    template <typename T, typename MemoryBlock = std::vector<T> >
    class general_matrix {
    public:
        // typedefs required for a std::container concept
        typedef T                       value_type;
        ///< The type T of the elements of the matrix
        typedef T&                      reference;
        ///< Reference to value_type
        typedef T const&                const_reference;
        ///< Const reference to value_type
        typedef std::size_t             size_type;
        ///< Unsigned integer type that represents the dimensions of the matrix
        typedef std::ptrdiff_t          difference_type;
        ///< Signed integer type to represent the distance of two elements in the memory

        // for compliance with an std::container one would also need
        // -operators == != < > <= >=
        // -size()
        // -typedefs iterator, const_iterator


        // typedefs for matrix specific iterators
        typedef strided_iterator<general_matrix,value_type>
            row_element_iterator;
        ///< Iterator to iterate through the elements of a row of the matrix
        typedef strided_iterator<const general_matrix,const value_type>
            const_row_element_iterator;
        ///< Const version of row_element_iterator
        typedef value_type*
            column_element_iterator;
        ///< Iterator to iterate through the elements of a columns of the matrix
        typedef value_type const*
            const_column_element_iterator;
        ///< Const version of column_element_iterator
        
        typedef matrix_element_iterator<general_matrix,value_type>
            element_iterator;
        ///< Iterator to iterate through all elements of the matrix (REALLY SLOW! USE row_-/column_iterators INSTEAD!)
        typedef matrix_element_iterator<const general_matrix,const value_type>
            const_element_iterator;
        ///< Const version of element_iterator (REALLY SLOW! USE row_-/column_iterators INSTEAD!)


        /**
          * The constructor
          * @param rows the number of rows
          * @param columns the number of columns
          * @param init_value all matrix elements will be initialized to this value.
          */
        general_matrix(size_type rows = 0, size_type columns = 0, T init_value = T() )
        : size1_(rows), size2_(columns), reserved_size1_(rows), values_(rows*columns, init_value)
        {
        }

        /**
          * The copy constructor
          *
          */
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
        
        /**
          * Non-throwing swap function
          * @param r general_matrix object which should be swapped with the general_matrix (this)
          */
        void swap(general_matrix & r)
        {
            std::swap(values_, r.values_);
            std::swap(size1_, r.size1_);
            std::swap(size2_, r.size2_);
            std::swap(reserved_size1_,r.reserved_size1_);
        }
        
        /**
          * Swaps two general_matrices
          */
        friend void swap(general_matrix & x, general_matrix & y)
        {
            x.swap(y);
        }

        /**
          * Assigns the matrix to matrix rhs
          */
        general_matrix& operator = (general_matrix rhs)
        {
            this->swap(rhs);
            return *this;
        }
        /**
          * Access the element in row i, column j
          * @param i 0<= i <= num_rows()
          * @param j 0<= j <= num_columns()
          * @return A mutable reference to the matrix element at position (i,j).
          */
        inline value_type& operator()(const size_type i, const size_type j)
        {
            assert((i < size1_) && (j < size2_));
            return values_[i+j*reserved_size1_];
        }
        
        /**
          * Access the element in row i, column j
          * @return A constant reference to the matrix element at position (i,j).
          */
        inline value_type const& operator()(const size_type i, const size_type j) const 
        {
            assert((i < size1_) && (j < size2_));
            return values_[i+j*reserved_size1_];
        }

        /**
          * Checks if a matrix is empty
          * @return true if the matrix is a 0x0 matrix, false otherwise.
          */
        inline const bool empty() const
        {
            return (size1_ == 0 || size2_ == 0);
        }

        /**
          * @return number of rows of the matrix
          */
        inline const size_type num_rows() const
        {
            return size1_;
        }

        /**
          * @return number of columns of the matrix
          */

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

        /**
          * @return The stride for moving to the next element along a row
          */
        inline const difference_type stride1() const
        {
            return 1;
        }
        
        /**
          * @return The stride for moving to the next element along a column
          */
        inline const difference_type stride2() const
        {
            return reserved_size1_;
        }

        /** Resize the matrix.
          * Resizes the matrix to size1 rows and size2 columns and initializes
          * the new elements to init_value. It also enlarges the MemoryBlock
          * if needed. If the new size for any dimension is
          * smaller only elements outside the new size will be deleted.
          * If the new size is larger for any dimension the new elements
          * will be initialized by the init_value.
          * All other elements will keep their value.  
          *
          * Exception behaviour:
          * As long as the assignment and copy operation of the T values don't throw an exception,
          * any exception will leave the matrix unchanged.
          * (Assuming the same behaviour of the underlying MemoryBlock. This is true for std::vector.)
          * @param size1 new number of rows
          * @param size2 new number of columns
          * @param init_value value to which the new elements will be initalized
          */
        void resize(size_type size1, size_type size2, T const& init_value = T())
        {

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
        /**
          * Reserves memory for anticipated enlargements of the matrix
          * @param size1 For how many rows should memory be reserved, value is ignored if it's smaller than the current number of rows
          * @param size2 For how many columns should memory be reserved, value is ignored if it's smaller than the current number of columns
          * @param init_value i
          */

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
        
        std::pair<element_iterator,element_iterator> elements()
        {
            return std::make_pair( element_iterator(this,0,0), element_iterator(this,0, num_columns()) );
        }
        
        std::pair<element_iterator,element_iterator> elements() const
        {
            return std::make_pair( const_element_iterator(this,0,0), const_element_iterator(this,0,num_columns() ) );
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
            for(size_type i=0; i< size1_; ++i)
                for(size_type j=0; j < size2_; ++j)
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
        
        void inplace_conjugate()
        {
            std::transform(this->elements().first, this->elements().second,
                           this->elements().first, functors::fconj());
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
    
    template<typename T, typename MemoryBlock>
    struct associated_diagonal_matrix<general_matrix<T, MemoryBlock> >
    {
        typedef diagonal_matrix<T> type;
    };
    
    
} // namespace blas

//
// Type promotion helper for mixed type matrix matrix and matrix vector operations
//
namespace blas {
template <typename T1, typename MemoryBlock1, typename T2, typename MemoryBlock2>
struct MultiplyReturnType
{
    typedef char one;
    typedef long unsigned int two;
    static one test(T1 t) {return one();}
    static two test(T2 t) {return two();}
    typedef typename boost::mpl::if_<typename boost::mpl::bool_<(sizeof(test(T1()*T2())) == sizeof(one))>,T1,T2>::type value_type;
    typedef typename boost::mpl::if_<typename boost::mpl::bool_<(sizeof(test(T1()*T2())) == sizeof(one))>,MemoryBlock1,MemoryBlock2>::type memoryblock_type;
};

// Specialize for U and T being the same type
template <typename T,typename MemoryBlock1, typename MemoryBlock2>
struct MultiplyReturnType<T,MemoryBlock1,T,MemoryBlock2>
{
    typedef T value_type;
    typedef MemoryBlock2 memoryblock_type;
};
}

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
            for(std::size_t k=0; k<lhs.num_columns(); ++k)
            {
                for(std::size_t j=0; j<rhs.num_columns(); ++j)
                {
                        result(i,j) += lhs(i,k) * rhs(k,j);
                }
            }
        }
        return result;
    } 
    
    template<typename T, typename MemoryBlock, typename T2, typename MemoryBlock2>
    const vector<typename MultiplyReturnType<T,MemoryBlock,T2,MemoryBlock2>::value_type,typename MultiplyReturnType<T,MemoryBlock,T2,MemoryBlock2>::memoryblock_type>
    matrix_vector_multiply(general_matrix<T,MemoryBlock> const& m, vector<T2,MemoryBlock2> const& v)
    {
        assert( m.num_columns() == v.size() );
        vector<
            typename MultiplyReturnType<T,MemoryBlock,T2,MemoryBlock2>::value_type,
            typename MultiplyReturnType<T,MemoryBlock,T2,MemoryBlock2>::memoryblock_type
            >
            result(m.num_rows());
        // Simple Matrix * Vector
        for(typename general_matrix<T,MemoryBlock>::size_type i = 0; i < m.num_rows(); ++i)
        {
            for(typename general_matrix<T,MemoryBlock>::size_type j=0; j <m.num_columns(); ++j)
            {
                result(i) += m(i,j) * v(j);
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

    template<typename T, typename MemoryBlock, typename T2, typename MemoryBlock2>
    const vector<typename MultiplyReturnType<T,MemoryBlock,T2,MemoryBlock2>::value_type, typename MultiplyReturnType<T,MemoryBlock,T2,MemoryBlock2>::memoryblock_type>
    operator * (general_matrix<T,MemoryBlock> const& m, vector<T2,MemoryBlock2> const& v)
    {
        return matrix_vector_multiply(m,v);
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
