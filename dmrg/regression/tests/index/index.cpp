#define BOOST_TEST_MODULE index
	
#include <numeric> 

#include <boost/geometry/geometries/adapted/boost_array.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/indexing.h"

using namespace index_detail;

typedef boost::mpl::list< U1, Zq<2> > SymmGroupList;

BOOST_AUTO_TEST_CASE_TEMPLATE(has_test, T, SymmGroupList){
    Index<T> phys;
    phys.insert(std::make_pair(0, 1));
    BOOST_CHECK_EQUAL(true,phys.has(0));
    BOOST_CHECK_EQUAL(false,phys.has(1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(size_of_block_test, T, SymmGroupList){
    Index<T> phys;
    phys.insert(std::make_pair(0, 2));
    BOOST_CHECK_EQUAL(2,phys.size_of_block(0));
    phys.insert(std::make_pair(1, 1));
    BOOST_CHECK_EQUAL(1,phys.size_of_block(1));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(positioni_test, T, SymmGroupList){
    Index<T> phys;
    phys.insert(std::make_pair(0, 5));
    phys.insert(std::make_pair(3, 2));
    phys.insert(std::make_pair(2, 3));
    phys.insert(std::make_pair(1, 4));
    phys.insert(std::make_pair(4, 1));
    BOOST_CHECK_EQUAL(1,phys.find(3));
    phys.insert(std::make_pair(3, 2));
    BOOST_CHECK_EQUAL(1,phys.find(3));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(position_pair_test, T, SymmGroupList){
    Index<T> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 2));
    phys.insert(std::make_pair(3, 1));
    BOOST_CHECK_EQUAL(2,phys.position(std::make_pair(2, 1)));
}

BOOST_AUTO_TEST_CASE(sort_test){
    Index<U1> phys;
    Index<U1> phys2;
    phys.insert(std::make_pair(2, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(0, 1));

    phys2.insert(std::make_pair(0, 1));
    phys2.insert(std::make_pair(1, 1));
    phys2.insert(std::make_pair(2, 1));

    BOOST_CHECK_EQUAL(phys,phys2);
}

BOOST_AUTO_TEST_CASE(operator_test){
    std::pair<U1::charge,std::size_t> p(std::make_pair(0,1)); 
    std::pair<U1::charge,std::size_t> p1(std::make_pair(1,1)); 

    BOOST_CHECK_EQUAL(true,p<p1);
    BOOST_CHECK_EQUAL(false,p1<p);

    Index<U1> phys;
    Index<U1> phys2;

    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(2, 1));

    phys2.insert(std::make_pair(0, 1));
    phys2.insert(std::make_pair(1, 1));
    phys2.insert(std::make_pair(2, 1));

    BOOST_CHECK_EQUAL(true,phys==phys2);
}

BOOST_AUTO_TEST_CASE(get_first_second_test){
    std::pair<U1::charge,std::size_t> p(std::make_pair(0,1)); 

    BOOST_CHECK_EQUAL(0,get_first<U1>(p));
    BOOST_CHECK_EQUAL(1,get_second<U1>(p));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(charges_test, T, SymmGroupList){
    Index<T> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    
    std::vector<typename T::charge> res(2);
    res[0] = 1; res[1] = 0;  
 
    std::vector<typename T::charge> vec = phys.charges();
   
    BOOST_CHECK_EQUAL(vec[0],res[0]);
    BOOST_CHECK_EQUAL(vec[1],res[1]);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sizes_test, T, SymmGroupList){
    Index<T> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 2));
    
    std::vector<std::size_t> res(2);
    res[0] = 2; res[1] = 1;  
 
    std::vector<std::size_t> vec = phys.sizes();
   
    BOOST_CHECK_EQUAL(vec[0],res[0]);
    BOOST_CHECK_EQUAL(vec[1],res[1]);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sum_of_size_est, T, SymmGroupList){
    Index<T> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 2));
    
    BOOST_CHECK_EQUAL(3,phys.sum_of_sizes());
}
