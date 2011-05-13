#include "utils/zout.hpp"
#include <vector>
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "p_dense_matrix/concept/matrix_interface.hpp"
#include "p_dense_matrix/concept/resizable_matrix_interface.hpp"

#define BOOST_TEST_MODULE p_memory_management
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <complex>
#include <numeric>

#define M_SIZE 6
using namespace blas;

typedef boost::mpl::list<double> test_types;
typedef ambient::dim2 dim;

struct caveats {
    caveats(){ ambient::init();     }
   ~caveats(){ ambient::finalize(); }
};

BOOST_GLOBAL_FIXTURE( caveats );

p_dense_matrix<double> foo(p_dense_matrix<double> a){

    a.resize(2,2);
    printf("Called foo!\n");
    p_dense_matrix<double> b = a;
    return a;
}


BOOST_AUTO_TEST_CASE_TEMPLATE( p_memory_management_test, T, test_types )
{
    ambient::layout >> dim(1,1), dim(1,1), dim(10,1);

    p_dense_matrix<double> d(4,4);
    p_dense_matrix<double> dc = foo(d);
    std::cout << d;
    //p_dense_matrix<double> e(0,0);
    //p_dense_matrix<double> f(0,0);
    //p_dense_matrix<double> dcc(d);
//    p_dense_matrix<double> dccc(dcc);

    //d.touch();
//    dcc.touch();

    //std::cout << d;
    //std::cout << dc;
    //dc.touch();
    //resize(e, 4, 4);
    //std::cout << e;
    //resize(f, 4, 4);
    //std::cout << f;

}

