#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.h"
#include "p_dense_matrix/concept/matrix_interface.hpp"
#include "p_dense_matrix/concept/resizable_matrix_interface.hpp"

#define BOOST_TEST_MODULE p_dense_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <complex>
#include <numeric>

#define M_SIZE 512
using namespace blas;

//
// List of types T for which the p_dense_matrix<T> is tested
// (long long unsigned int causes problems in boost::iterator facade)
typedef boost::mpl::list<double> test_types;
typedef ambient::dim3 dim3;

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

struct AmbientConfig {
    AmbientConfig()   { ambient::init();     }
    ~AmbientConfig()  { ambient::finalize(); }
};


//
// List of type pairs <T,U> for which the mixed type matrix vector multiplication is tested.
//
typedef boost::mpl::list<type_pairs::IntDouble, type_pairs::DoubleInt> test_type_pairs;

BOOST_GLOBAL_FIXTURE( AmbientConfig );

/*BOOST_AUTO_TEST_CASE_TEMPLATE( summ_operation_test, T, test_types )
{
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);
    p_dense_matrix<T> d(M_SIZE,M_SIZE);

    c = a + b;
    c(5,5) = 13.0;
    c(6,5) = 14.0;
    a.remove_rows(128,128);
    b.remove_rows(128,128);
    c.remove_rows(128,128);
    c = a + b;

    a.resize(640,512);
    b.resize(640,512);
    c.resize(640,512);
    c = a + b;
    ambient::playout();
}*/

BOOST_AUTO_TEST_CASE_TEMPLATE( sql_test, T, test_types )
{
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);

    c = a + b;
    c = a + b;
    ambient::playout();
}

