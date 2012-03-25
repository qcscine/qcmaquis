#define BOOST_TEST_MODULE vli_polynomial
#include <boost/test/unit_test.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/static_assert.hpp>

#include "vli/polynomial/monomial.hpp"
#include "vli/polynomial/polynomial_cpu.h"
#include "vli/vli_cpu.h"

#include "regression/vli_test.hpp"

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::vli_cpu;
using vli::monomial;
using vli::polynomial_cpu;

template <typename Vli>
struct polynomial_from_vli_cpu
{
    typedef polynomial_cpu<Vli, 10> type;
};

typedef boost::mpl::transform<
      vli::test::vli_cpu_type_list
    , polynomial_from_vli_cpu<boost::mpl::_1>
    >::type polynomial_types;
        

BOOST_AUTO_TEST_CASE_TEMPLATE( construction_and_coeff_assignment, Poly, polynomial_types )
{
    Poly pa;
    typedef typename Poly::exponent_type exponent_type;
    typename Poly::value_type a[Poly::max_order*Poly::max_order];

    for(exponent_type i=0; i < Poly::max_order; ++i)
        for(exponent_type j=0; j < Poly::max_order; ++j)
            BOOST_CHECK_EQUAL(pa(i,j), typename Poly::value_type(0));

    for(exponent_type i=0; i < Poly::max_order; ++i)
    {
        for(exponent_type j=0; j < Poly::max_order; ++j)
        {
            fill_random(a[i*Poly::max_order+j]);
            pa(i,j) = a[i*Poly::max_order+j];
        }
    }

    for(exponent_type i=0; i < Poly::max_order; ++i)
        for(exponent_type j=0; j < Poly::max_order; ++j)
            BOOST_CHECK_EQUAL(pa(i,j),a[i*Poly::max_order+j]);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( copy_constructor_and_equal, Poly, polynomial_types )
{
    BOOST_STATIC_ASSERT(Poly::max_order >= 2);
    Poly pa;
    fill_poly_random(pa);

    Poly pb(pa);

    for(typename Poly::value_type::size_type i=0; i < Poly::max_order; ++i)
        for(typename Poly::value_type::size_type j=0; j < Poly::max_order; ++j)
            BOOST_CHECK_EQUAL(pb(i,j),pa(i,j));

    BOOST_CHECK_EQUAL(pa == pb, true);

    pb(Poly::max_order-2, Poly::max_order-2) += 1;

    BOOST_CHECK_EQUAL(pa == pb, false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignement_operator, Poly, polynomial_types )
{
    BOOST_STATIC_ASSERT(Poly::max_order >= 2);
    Poly pa;
    fill_poly_random(pa);

    Poly pb;
    fill_poly_random(pb);

    BOOST_CHECK_EQUAL( pa == pb, false);

    pb = pa;

    BOOST_CHECK_EQUAL(pb, pa);
    
    // Check whether pa changes when pb is changed.
    pb(1,1) = typename Poly::value_type(857);

    BOOST_CHECK_EQUAL(pb == pa, false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign, Poly, polynomial_types )
{
    Poly pa;
    Poly pb;
    fill_poly_random(pb);

    pa += pb;

    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign_zeroth_order, Poly, polynomial_types )
{
    Poly pa;
    fill_poly_random(pa);

    Poly pb(pa);

    int a(3);
    pa += a;

    BOOST_CHECK_EQUAL(pa(0,0), pb(0,0)+a);
    
    pa(0,0) = pb(0,0);
    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign_monomial, Poly, polynomial_types )
{
    Poly pa;
    fill_poly_random(pa);
    typename Poly::value_type va;
    fill_random(va);

    Poly pb(pa);

    monomial<typename Poly::value_type> ma(1,0);
    ma *= va;
    pa += ma;
    pb(1,0) += va;
    BOOST_CHECK_EQUAL(pa,pb);

    typename Poly::value_type vb;
    fill_random(vb);
    monomial<typename Poly::value_type> mb(0,1);
    mb *= vb;
    pa += mb;
    pb(0,1) += vb;
    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign, Poly, polynomial_types )
{
    Poly pa;
    Poly pb;

    fill_poly_random(pb);

    pa -= pb;

    for(typename Poly::value_type::size_type i=0; i < Poly::max_order; ++i)
        for(typename Poly::value_type::size_type j=0; j < Poly::max_order; ++j)
            BOOST_CHECK_EQUAL(pa(i,j), -pb(i,j));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign_zeroth_order, Poly, polynomial_types )
{
    Poly pa;
    fill_poly_random(pa);

    Poly pb(pa);

    int a=9;
    pa -= a;

    BOOST_CHECK_EQUAL(pa(0,0), pb(0,0)-a);
    
    pa(0,0) = pb(0,0);
    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign_monomial, Poly, polynomial_types )
{
    Poly pa;
    fill_poly_random(pa);
    typename Poly::value_type va;
    fill_random(va);

    Poly pb(pa);

    monomial<typename Poly::value_type> ma(1,0);
    ma *= va;
    pa -= ma;
    pb(1,0) -= va;
    BOOST_CHECK_EQUAL(pa,pb);

    
    typename Poly::value_type vb;
    fill_random(vb);
    monomial<typename Poly::value_type> mb(0,1);
    mb *= vb;
    pa -= mb;
    pb(0,1) -= vb;
    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_assign_constant, Poly, polynomial_types )
{
    Poly pa;
    fill_poly_random(pa);

    Poly pb(pa);
    typename Poly::value_type a;
    fill_random(a);

    pa*=a;

    for(typename Poly::value_type::size_type i=0; i < Poly::max_order; ++i)
        for(typename Poly::value_type::size_type j=0; j < Poly::max_order; ++j)
            BOOST_CHECK_EQUAL(pa(i,j),pb(i,j)*a);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_monomial, Poly, polynomial_types )
{
    Poly pa;
    fill_poly_random(pa);

    typename Poly::value_type vli_zero(0);
    typename Poly::value_type va(3);
    monomial<typename Poly::value_type> ma(1,0);
    ma *= va;

    Poly pb = pa * ma;

    BOOST_CHECK_EQUAL(pb(0,0),vli_zero);
    BOOST_CHECK_EQUAL(pb(0,1),vli_zero);
    for(typename Poly::value_type::size_type i=0; i < Poly::max_order-1; ++i)
        for(typename Poly::value_type::size_type j=0; j < Poly::max_order; ++j)
            BOOST_CHECK_EQUAL(pb(i+1,j), pa(i,j)*va);

    typename Poly::value_type vb(5);
    monomial<typename Poly::value_type> mb(0,1);
    mb *= vb;
    
    pb = pb * mb;
    BOOST_CHECK_EQUAL(pb(0,0),vli_zero);
    BOOST_CHECK_EQUAL(pb(0,1),vli_zero);
    BOOST_CHECK_EQUAL(pb(1,0),vli_zero);
    for(typename Poly::value_type::size_type i=0; i < Poly::max_order-1; ++i)
        for(typename Poly::value_type::size_type j=0; j < Poly::max_order-1; ++j)
            BOOST_CHECK_EQUAL(pa(i,j)*va*vb, pb(i+1,j+1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies, Poly, polynomial_types )
{    
    typedef vli::vli_cpu<typename  Poly::value_type::size_type, 2*(Poly::value_type::size) > vli_result_type_cpu;
    typedef vli::vli_cpu<typename  Poly::value_type::size_type, (Poly::value_type::size) > vli_type_cpu;    
    typedef vli::polynomial_cpu<vli_result_type_cpu, 2*(Poly::max_order) > polynomial_result_type_cpu;
    typedef vli::polynomial_cpu<vli_type_cpu, (Poly::max_order) > polynomial_type_cpu;
    
    polynomial_result_type_cpu pc;
    polynomial_type_cpu pa;
    fill_poly_random(pa);

    polynomial_type_cpu pb;
    fill_poly_random(pb);
   // poly_multiply(pc,pa,pb);
    pc = pa*pb;

    polynomial_result_type_cpu result;
    for(typename Poly::value_type::size_type i=0; i < Poly::max_order; ++i)
        for(typename Poly::value_type::size_type j=0; j < Poly::max_order; ++j)
            for(typename Poly::value_type::size_type k=0; k < (Poly::max_order ); ++k)
                for(typename Poly::value_type::size_type l=0; l < (Poly::max_order); ++l)
                    muladd(result(i+k,j+l),pa(i,j),pb(k,l));

    BOOST_CHECK_EQUAL(pc,result);    

 }

