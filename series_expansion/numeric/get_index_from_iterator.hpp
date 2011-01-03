#ifndef HP2C__GET_INDEX_FROM_ITERATOR_HPP
#define HP2C__GET_INDEX_FROM_ITERATOR_HPP

#include <cassert>
#include <vector>

namespace series_expansion
{

    /*
template <typename Vector>
inline typename Vector::size_type get_index_from_iterator(Vector const& v, typename Vector::const_iterator it)
{
    std::runtime_error("get_index_from_iterator() is not implemeneted for the used Vector type.");
    return typename Vector::size_type(0);
}
*/

#ifdef HP2C__SIMPLE_SPARSE_VECTOR_HPP
template <typename T>
inline typename simple_sparse_vector<T>::size_type get_index_from_iterator(simple_sparse_vector<T> const& v, typename simple_sparse_vector<T>::const_iterator it)
{
    //assert( v.begin() <= it);
    //assert( it < v.end() );
    return it.get_index();
}
#endif //HP2C__SIMPLE_SPARSE_VECTOR_HPP

// Specialization for std::vector
template <typename T>
inline typename std::vector<T>::size_type get_index_from_iterator(std::vector<T> const& v, typename std::vector<T>::const_iterator it)
{
    assert( v.begin() <= it );
    assert( it < v.end() );
    return std::distance( v.begin(), it );
}


}

#endif //HP2C__GET_INDEX_FROM_ITERATOR_HPP
