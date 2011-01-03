#ifndef HP2C__SIMPLE_SPARSE_VECTOR_ITERATOR_HPP
#define HP2C__SIMPLE_SPARSE_VECTOR_ITERATOR_HPP

#include <boost/iterator/iterator_facade.hpp>
#include <boost/static_assert.hpp>
#include <cassert>

namespace series_expansion
{

template <typename T>
    class simple_sparse_vector;

template <typename T, typename MapIterator, typename Index>
class simple_sparse_vector_iterator : public boost::iterator_facade<
                                simple_sparse_vector_iterator<T,MapIterator,Index>,
                                T,
                                boost::bidirectional_traversal_tag,
                                T&
                                >
{
    public:
        typedef T value_type;

        simple_sparse_vector_iterator(MapIterator m_it)
            : m_it(m_it)
        {
            // The value_type of the iterator must be the value_type of the matrix or const Matrix::value_type
//            BOOST_STATIC_ASSERT( (boost::is_same<typename Matrix::value_type, T>::value
//                                 || boost::is_same<const typename Matrix::value_type,T>::value) );
        }

        template<typename T2, typename MapIterator2>
        simple_sparse_vector_iterator(simple_sparse_vector_iterator<T2,MapIterator2,Index> const& r)
            : m_it(r.m_it)
        {
            BOOST_STATIC_ASSERT( (boost::is_same<
                        typename boost::add_const<T>::type,
                        typename boost::add_const<T2>::type
                        >::value ));
        }

        Index get_index() const
            { return m_it->first; }
    private:
        friend class boost::iterator_core_access;
       
        template <typename T2, typename MapIterator2, typename Index2>
        friend class simple_sparse_vector_iterator;

        value_type& dereference() const
            { return m_it->second; }

        template <typename T2, typename MapIterator2>
        bool equal(simple_sparse_vector_iterator<T2,MapIterator2,Index> const& y) const
        {
            BOOST_STATIC_ASSERT( (boost::is_same<
                        typename boost::add_const<T>::type,
                        typename boost::add_const<T2>::type
                        >::value ));
            if(m_it == y.m_it)
                return true;
            else
                return false;
        }

        void increment()
        {
            ++m_it;
        }

        void decrement()
        {
            --m_it;
        }

//        template <typename T2, typename MapIterator2>
//        std::ptrdiff_t distance_to(sparse_vector_iterator<T2,MapIterator2,Index> const& y) const
//        {
//            std::cerr<<"WARNING distance_to for the sparse_vector_iterator is very slow"<<std::endl;
//            return std::distance(y.m_it,m_it);
//        }

        // The actual element
        MapIterator m_it;
};
}
#endif //HP2C__SIMPLE_SPARSE_VECTOR_ITERATOR_HPP
