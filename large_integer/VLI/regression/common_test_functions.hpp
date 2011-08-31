#ifndef COMMON_TEST_FUNCTIONS_HPP
#define COMMON_TEST_FUNCTIONS_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <vli_cpu/vli_number_traits.hpp>

namespace vli
{
namespace test
{
boost::mt11213b rng;

template <typename Vli>
typename Vli::value_type rnd_digit()
{
    static boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);
    return rnd(rng);
}

template <typename Vli>
int rnd_valid_int()
{
    static boost::uniform_int<int> rnd(0,std::abs(static_cast<int>(max_int_value<Vli>::value)));
    return rnd(rng);
}

template <typename Vli>
void fill_random(Vli& v)
{
    typename Vli::value_type midle = Vli::size/2-1;

    for(typename Vli::size_type i=0; i < midle; ++i)
        v[i] = rnd_digit<Vli>();
}

template <typename Vli>
void fill_random(Vli& v, typename Vli::size_type size)
{
    assert(size <= Vli::size);
    for(typename Vli::size_type i=0; i < (size); ++i)
        v[i] = rnd_digit<Vli>();
}

template <typename Polynomial>
void fill_poly_random(Polynomial& p)
{
    for(typename Polynomial::size_type i=0; i < Polynomial::max_order; ++i)
        for(typename Polynomial::size_type j=0; j < Polynomial::max_order; ++j)
            fill_random(p(i,j));
}

template <typename Vector>
void fill_vector_random(Vector& v)
{
    for(typename Vector::size_type i=0; i < v.size(); ++i)
        fill_poly_random(v[i]);
}



} //namespace test
} //namespace vli

#endif //COMMON_TEST_FUNCTIONS_HPP
