#define BOOST_TEST_MODULE flexible_polynomial
#include <boost/test/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/size.hpp>
#include <limits>
#include <fstream>
#include <iostream>

#include <vli/polynomial/polynomial.hpp>

using vli::var;
typedef boost::mpl::vector<
      vli::polynomial<int,vli::max_order_each<4>,     var<'x'> >
    , vli::polynomial<int,vli::max_order_each<4>,     var<'x'>, var<'z'> >
    , vli::polynomial<int,vli::max_order_each<4>,     var<'x'>, var<'y'>, var<'z'>, var<'a'> >
    , vli::polynomial<int,vli::max_order_combined<4>, var<'y'> >
    , vli::polynomial<int,vli::max_order_combined<4>, var<'x'>, var<'y'> >
    , vli::polynomial<int,vli::max_order_combined<4>, var<'x'>, var<'y'>, var<'z'> >
    , vli::polynomial<int,vli::max_order_combined<4>, var<'x'>, var<'y'>, var<'z'>, var<'a'> >
    > polynomial_types;

struct reference_data_storage {
    reference_data_storage() {

        // TODO think of something more solid
        std::fstream ref_file("polynomial_tests_reference.txt");

        if(!ref_file.good())
        {
            ref_file.close();
            ref_file.open("regression/polynomial_tests_reference.txt");
        }
        if(!ref_file.good())
            throw std::runtime_error("Reference data file not found!");

        while(ref_file.good()) {
            std::string line;
            std::getline(ref_file,line);
            if(line[0] != '#')
                ref_strings.push_back(line);
        }
        ref_file.close();
        std::cout<<ref_strings.size()<<" reference data sets loaded."<<std::endl;
    }
    template <class TestType>
    std::string get_result(unsigned int test_index, TestType type) const {
        typedef typename boost::mpl::find<polynomial_types,TestType>::type mpl_iter;
        unsigned int type_index = mpl_iter::pos::value;
        std::size_t index = test_index * boost::mpl::size<polynomial_types>::value + type_index;
        if(index >= ref_strings.size() )
            throw std::runtime_error("Reference data not found.");
        return ref_strings[index];
    }
 private:
   std::vector<std::string> ref_strings;
};

static const reference_data_storage reference_data;


template <typename Polynomial>
int fill_polynomial(Polynomial& p,int offset=0) {
    for(typename Polynomial::iterator it=p.begin(); it != p.end(); ++it)
        *it = offset++;
    return offset;
}


BOOST_AUTO_TEST_CASE_TEMPLATE ( constructor, Poly, polynomial_types ) {
    typedef typename Poly::value_type value_type;

    Poly a;
    for(typename Poly::iterator it=a.begin(); it != a.end(); ++it)
        BOOST_CHECK_EQUAL(*it,value_type(0));
    
    value_type x = std::numeric_limits<value_type>::max();
    Poly b(x);
    typename Poly::iterator it = b.begin();
    BOOST_CHECK_EQUAL(*it,x);
    ++it;
    for(; it != b.end(); ++it)
        BOOST_CHECK_EQUAL(*it,value_type(0));
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( copy_constructor, Poly, polynomial_types ) {
    Poly a;
    fill_polynomial(a);

    Poly b(a);

    typename Poly::iterator it(a.begin()), it2(b.begin());
    while(it != a.end())
        BOOST_CHECK_EQUAL(*it++,*it2++);
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( converting_copy_constructor, Poly, polynomial_types ) {
    typedef vli::polynomial<
          double
        , typename Poly::max_order
        , typename vli::variable<Poly,0>::type
        , typename vli::variable<Poly,1>::type
        , typename vli::variable<Poly,2>::type
        , typename vli::variable<Poly,3>::type
        > polynomial_double_type;
    Poly a;
    fill_polynomial(a);
    polynomial_double_type b(a);

    typename Poly::iterator                     it(a.begin());
    typename polynomial_double_type::iterator   it2(b.begin());
    while(it != a.end())
        BOOST_CHECK_EQUAL(static_cast<double>(*it++),*it2++);
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( equal_constructor, Poly, polynomial_types ) {
    Poly a;
    fill_polynomial(a);
    Poly b(a);

    BOOST_CHECK_EQUAL(a == b, true);

    // Change one element
    assert(std::distance(a.begin(),a.end()) > 3);
    *(a.end()-3) += 1;

    BOOST_CHECK_EQUAL(a == b, false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignement_operator, Poly, polynomial_types ) {
//    BOOST_STATIC_ASSERT(Poly::max_order >= 2);
    Poly pa;
    Poly pb;
    int x = fill_polynomial(pa);
    fill_polynomial(pb,x);

    BOOST_CHECK_EQUAL( pa == pb, false);

    pb = pa;

    BOOST_CHECK_EQUAL(pb, pa);
    
    // Check whether pa changes when pb is changed.
    *(pb.end()-3) += 1;

    BOOST_CHECK_EQUAL(pb == pa, false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( print, Poly, polynomial_types ) {
    Poly a;
    fill_polynomial(a);
    std::stringstream a_str;
    a_str<<a;
    BOOST_CHECK_EQUAL(a_str.str(),(reference_data.get_result(0,a)));
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( swap_test, Poly, polynomial_types ) {
    Poly pa;
    Poly pb;

    int x = fill_polynomial(pa,0);
    fill_polynomial(pb,x);
    
    Poly pc(pa);
    Poly pd(pb);

    BOOST_CHECK_EQUAL(pa,pc);
    BOOST_CHECK_EQUAL(pb,pd);

    swap(pa,pb);

    BOOST_CHECK_EQUAL(pb,pc);
    BOOST_CHECK_EQUAL(pa,pd);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign, Poly, polynomial_types) {
    Poly pa;
    Poly pb;
    fill_polynomial(pb);


    pa += pb;

    BOOST_CHECK_EQUAL(pa,pb);

    fill_polynomial(pb,63243);

    for(typename Poly::iterator it=pa.begin(); it != pa.end(); ++it)
        *it = 1;

    pa += pb;

    typename Poly::iterator it(pa.begin()), it2(pb.begin());
    while(it != pa.end())
        BOOST_CHECK_EQUAL(*it++,(*it2++)+1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign_zeroth_order, Poly, polynomial_types) {
    Poly pa;
    fill_polynomial(pa);

    Poly pb(pa);

    int a(3);
    pa += a;

    BOOST_CHECK_EQUAL(pa(0), pb(0)+a);
    
    pa(0) = pb(0);
    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign_monomial, Poly, polynomial_types) {
    // TODO improve by checking the 2nd, 3rd or 4th variable too
    Poly pa;
    fill_polynomial(pa);
    Poly pb(pa);
    typename Poly::value_type va = 5;

    vli::monomial<typename Poly::value_type, typename vli::variable<Poly,0>::type> ma(2);
    ma *= va;
    pa += ma;
    pb(2) += va;
    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign, Poly, polynomial_types ) {
    Poly pa;
    Poly pb;
    fill_polynomial(pb);

    pa -= pb;

    typename Poly::iterator it(pa.begin()), it2(pb.begin());
    while( it != pa.end() )
        BOOST_CHECK_EQUAL(*it++,-*it2++);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign_zeroth_order, Poly, polynomial_types ) {
    Poly pa;
    fill_polynomial(pa);
    Poly pb(pa);

    int a=9;
    pa -= a;

    BOOST_CHECK_EQUAL(pa(0), pb(0)-a);
    
    pa(0) = pb(0);
    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign_monomial, Poly, polynomial_types ) {
    // TODO improve by checking the 2nd, 3rd or 4th variable too
    Poly pa;
    fill_polynomial(pa);
    Poly pb(pa);
    typename Poly::value_type va = 5;

    vli::monomial<typename Poly::value_type, typename vli::variable<Poly,0>::type> ma(2);
    ma *= va;
    pa -= ma;
    pb(2) -= va;
    BOOST_CHECK_EQUAL(pa,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_assign_constant, Poly, polynomial_types) {
    Poly pa;
    fill_polynomial(pa);
    Poly pb(pa);

    typename Poly::value_type a = 42;
    pa *= a;

    typename Poly::iterator it(pa.begin()), it2(pb.begin());
    while( it != pa.end() )
        BOOST_CHECK_EQUAL(*it++,(*it2++)*a);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_assign_int, Poly, polynomial_types) {
    Poly pa;
    fill_polynomial(pa);
    Poly pb(pa);

    int a = 42;
    pa *= a;

    typename Poly::iterator it(pa.begin()), it2(pb.begin());
    while( it != pa.end() )
        BOOST_CHECK_EQUAL(*it++,(*it2++)*a);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( negate, Poly, polynomial_types ) {
    Poly pa;
    fill_polynomial(pa);
    Poly pb(pa);

    negate_inplace(pa);

    typename Poly::iterator it(pa.begin()), it2(pb.begin());
    while( it != pa.end() )
        BOOST_CHECK_EQUAL(*it++,-*it2++);
    
    it = pa.begin(); it2 = pb.begin();
    while( it != pa.end() )
        BOOST_CHECK_EQUAL(*it++,typename Poly::value_type(0)-(*it2++));

    Poly pc = -pa;

   BOOST_CHECK_EQUAL(pc,pb);
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( monomial_multiply, Poly, polynomial_types ) {
    // TODO improve by checking the 2nd, 3rd or 4th variable too
    Poly pa;
    fill_polynomial(pa);
    vli::monomial<typename Poly::value_type, typename vli::variable<Poly,0>::type> m(2);
    m *= 10;
    pa *= m;

    std::stringstream pa_str;
    pa_str<<pa;
    BOOST_CHECK_EQUAL(pa_str.str(),(reference_data.get_result(1,pa)));
}

BOOST_AUTO_TEST_CASE_TEMPLATE ( polynomial_multiply, Poly, polynomial_types ) {
    Poly pa;
    Poly pb;
    int x = fill_polynomial(pa,0);
    fill_polynomial(pb,x);
    
    std::stringstream papb_str;
    papb_str<<(pa*pb);
    BOOST_CHECK_EQUAL(papb_str.str(),(reference_data.get_result(2,pa)));
}

