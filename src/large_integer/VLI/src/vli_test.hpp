#ifndef COMMON_TEST_FUNCTIONS_HPP
#define COMMON_TEST_FUNCTIONS_HPP

#include <boost/mpl/list.hpp>
#include <boost/mpl/pop_front.hpp>

#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

#include "vli/vli_config.h"
#include "vli/vli_traits.hpp"


namespace vli
{
namespace test
{


//
// Generate a list of all vli types for which the library was compiled
// and run tests for them.
// (VLI_COMPILE_BASEINT_SIZE_PAIRS_SEQ is defined by CMake
//  in vli_utils/vli_config.h)
//

#define VLI_CPU_APPEND_TEST_TYPE(r,data,BASEINT_SIZE_PAIR) \
        , vli_cpu< BOOST_PP_TUPLE_ELEM(2,0,BASEINT_SIZE_PAIR) , BOOST_PP_TUPLE_ELEM(2,1,BASEINT_SIZE_PAIR) >

typedef boost::mpl::pop_front<
    boost::mpl::list<
        boost::mpl::void_
BOOST_PP_SEQ_FOR_EACH(VLI_CPU_APPEND_TEST_TYPE, _, VLI_COMPILE_BASEINT_SIZE_PAIRS_SEQ)
    >
>::type vli_cpu_type_list;



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
    for(typename Vli::size_type i=0; i < Vli::size; ++i)
        v[i] = rnd_digit<Vli>();
}

template <typename Vli>
void fill_random(Vli& v, typename Vli::size_type size)
{
    assert(size <= Vli::size);
    for(typename Vli::size_type i=0; i < size; ++i)
        v[i] = rnd_digit<Vli>();
}

template <typename Polynomial>
void fill_poly_random(Polynomial& p)
{
    for(typename Polynomial::exponent_type i=0; i < Polynomial::max_order; ++i)
        for(typename Polynomial::exponent_type j=0; j < Polynomial::max_order; ++j)
            fill_random(p(i,j),1);
}

template <typename Polynomial>
void fill_poly_random(Polynomial& p, typename Polynomial::exponent_type size)
{
    for(typename Polynomial::exponent_type i=0; i < Polynomial::max_order; ++i)
        for(typename Polynomial::exponent_type j=0; j < Polynomial::max_order; ++j)
            fill_random(p(i,j),size);
}

template <typename Vector>
void fill_vector_random(Vector& v)
{
    for(typename Vector::size_type i=0; i < v.size(); ++i)
        fill_poly_random(v[i]);
}

template <typename Vector>
void fill_vector_random(Vector& v, typename Vector::value_type::exponent_type size)
{
    for(typename Vector::size_type i=0; i < v.size(); ++i)
        fill_poly_random(v[i], size);
}



} //namespace test
} //namespace vli

#endif //COMMON_TEST_FUNCTIONS_HPP
