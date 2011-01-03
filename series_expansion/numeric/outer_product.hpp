#ifndef HP2C__OUTER_PRODUCT_HPP
#define HP2C__OUTER_PRODUCT_HPP

#include <cassert>


namespace blas
{

// TODO
// - until now it is a single special case of a generalized outer product
// - return type deduction
template <typename Vector>
const dense_matrix<typename Vector::value_type> outer_product(blas::vector<Vector> const& v1, blas::vector<Vector> const& v2)
{
    assert( v1.size() == v2.size() );
    dense_matrix<typename Vector::value_type> result(v1.size(),v2.size());
    for(unsigned int i = 0; i<v1.size(); ++i)
    {
        for(unsigned int j = 0; j<v2.size(); ++j)
        {
            result(i,j) = inner_product( v1[i], v2[j] );
        }
    }
    return result; 
}

}

#endif //HP2C__OUTER_PRODUCT_HPP
