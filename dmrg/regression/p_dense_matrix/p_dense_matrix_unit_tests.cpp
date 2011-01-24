#define BOOST_TEST_MODULE p_dense_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <complex>
#include <numeric>

#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.h"
#include "p_dense_matrix/concept/matrix_interface.hpp"
#include "p_dense_matrix/concept/resizable_matrix_interface.hpp"

using namespace blas;
//
// List of types T for which the p_dense_matrix<T> is tested
//
typedef boost::mpl::list<double, int> test_types;
// long long unsigned int causes problems in boost::iterator facade

namespace type_pairs
{

struct IntDouble
{
    typedef int first_type;
    typedef double second_type;
};

struct DoubleInt
{
    typedef double first_type;
    typedef int second_type;
};
};

//
// List of type pairs <T,U> for which the mixed type matrix vector multiplication is tested.
//
typedef boost::mpl::list<type_pairs::IntDouble, type_pairs::DoubleInt> test_type_pairs;

template <typename OutputIterator, typename T>
T fill_range_with_numbers(OutputIterator begin, OutputIterator end, T iota)
{
    // Unfortunately we can't use the postincrement operator, due to std:complex<>
    // -> so we have to emulate it's behaviour...
    std::transform(begin,end,begin,boost::lambda::_1 = (boost::lambda::var(iota)+=T(1))-T(1));
    return iota;
}

template <typename T>
T fill_matrix_with_numbers(p_dense_matrix<T>& a)
{
    T iota(0);
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename p_dense_matrix<T>::row_element_iterator, typename p_dense_matrix<T>::row_element_iterator> range(row(a,i));
        iota += fill_range_with_numbers(range.first,range.second,T(i));
    }
    return iota;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors_test, T, test_types )
{
    p_dense_matrix<T> a;
    BOOST_CHECK_EQUAL(num_rows(a), 0 );

    BOOST_CHECK_EQUAL(num_columns(a), 0 );

    p_dense_matrix<T> b(10,10);
    BOOST_CHECK_EQUAL(num_rows(b), 10 );
    BOOST_CHECK_EQUAL(num_columns(b), 10 );
    for(unsigned int i=0; i<10; ++i)
        for(unsigned int j=0; j<10; ++j)
            BOOST_CHECK_EQUAL(b(i,j), T());

    p_dense_matrix<T> c(15,5,5);
    BOOST_CHECK_EQUAL(num_rows(c), 15 );
    BOOST_CHECK_EQUAL(num_columns(c), 5 );
    for(unsigned int i=0; i<15; ++i)
        for(unsigned int j=0; j<5; ++j)
            BOOST_CHECK_EQUAL(c(i,j), T(5));
}



