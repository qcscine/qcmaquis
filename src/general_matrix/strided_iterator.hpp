#ifndef ALPS_STRIDED_ITERATOR
#define ALPS_STRIDED_ITERATOR

#include <boost/iterator/iterator_facade.hpp>
#include <boost/static_assert.hpp>

template <typename Matrix, typename T>
class strided_iterator : public boost::iterator_facade<
                                strided_iterator<Matrix,T>,
                                T,
                                boost::random_access_traversal_tag,
                                T&,
                                typename Matrix::difference_type
                                >
{
    public:
        typedef T value_type;

        strided_iterator(value_type* ptr, typename Matrix::difference_type stride)
            : ptr(ptr), stride(stride)
        {
            // The value_type of the iterator must be the value_type of the matrix or const Matrix::value_type
            BOOST_STATIC_ASSERT( (boost::is_same<typename Matrix::value_type, T>::value
                                 || boost::is_same<const typename Matrix::value_type,T>::value) );
        }

        template<typename Matrix2, typename U>
        strided_iterator(strided_iterator<Matrix2,U> const& r)
            : ptr(r.ptr), stride(r.stride)
            {}

    private:
        friend class boost::iterator_core_access;
        template <typename,typename> friend class strided_iterator;

        value_type& dereference() const
        { return *ptr; }

        // iterators are equal if they point to the same column of the same matrix
        // WARNING: since the row position is not compared
        // two iterators can be equal although they point to different elements
        template <typename U>
        bool equal(strided_iterator<Matrix,U> const& y) const
        {
            if(ptr == y.ptr)
                return true;
            else
                return false;
        }
        void increment()
        {
            ptr+=stride;
        }
        void decrement()
        {
            ptr-=stride;
        }
        void advance(typename Matrix::difference_type n)
        {
            ptr += n*stride;
        }

        template <typename U>
        typename Matrix::difference_type distance_to(strided_iterator<Matrix,U> const& z) const
        {
            assert( (ptr - z.ptr) % stride == 0 );
            return (ptr - z.ptr)/stride;
        }

        value_type* ptr;
        typename Matrix::difference_type stride; 
};


#endif //ALPS_STRIDED_ITERATOR
