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
    AmbientConfig()   { ambient::init(); }
    ~AmbientConfig()  { ambient::finalize(); }
};


//
// List of type pairs <T,U> for which the mixed type matrix vector multiplication is tested.
//
typedef boost::mpl::list<type_pairs::IntDouble, type_pairs::DoubleInt> test_type_pairs;

BOOST_GLOBAL_FIXTURE( AmbientConfig );

BOOST_AUTO_TEST_CASE_TEMPLATE( summ_operation_test, T, test_types )
{
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);
    p_dense_matrix<T> d(M_SIZE,M_SIZE);

    ambient::push(ambient::mem_bound_l_kernel, ambient::null_c_kernel, a, b, c);
    c = a + b;
    c(5,5) = 13.0;
    c(6,5) = 14.0;
    c.remove_rows(5,1);
    c.remove_cols(5,1);
    for(int i=0; i < 20; i++) if(c(i,5) >= 1.0) 
    printf("The element of c(%d,%d) is %.2f\n", i, 5, c(i,5));

////////////////////
    printf("\n\n\n\n\n\n\n\n");
    int*const *num = new int*((int*)malloc(sizeof(int)*10));
    for(int i = 0; i < 10; i++) (*num)[i] = i+13;
    ambient::push(ambient::single_integer_l_kernel, ambient::single_integer_c_kernel, *num);
    ambient::playout();


}

