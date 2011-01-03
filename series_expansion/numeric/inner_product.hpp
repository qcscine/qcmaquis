#ifndef HP2C__INNER_PRODUCT_HPP
#define HP2C__INNER_PRODUCT_HPP

#include <cassert>

namespace series_expansion
{

#ifdef HP2C__SIMPLE_SPARSE_VECTOR_HPP
template <typename T>
T inner_product(simple_sparse_vector<T> const& v1, simple_sparse_vector<T> const& v2)
{
    //assert(v1.size() == v2.size());
    return v1*v2;
}
#endif //HP2C__SIMPLE_SPARSE_VECTOR_HPP

template <typename T>
T inner_product(std::vector<T> const& v1, std::vector<T> const& v2)
{
    assert(v1.size() == v2.size());
    return std::inner_product(v1.begin(),v1.end(),v2.begin(),T(0) );
}

}

#endif //HP2C__INNER_PRODUCT_HPP
