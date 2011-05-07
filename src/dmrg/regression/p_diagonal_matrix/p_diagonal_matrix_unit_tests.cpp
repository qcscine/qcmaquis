#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_diagonal_matrix.h"

#include "p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "p_dense_matrix/concept/matrix_interface.hpp"
#include "p_dense_matrix/concept/resizable_matrix_interface.hpp"

#define BOOST_TEST_MODULE p_dense_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <complex>
#include <numeric>

#define M_SIZE 16
using namespace blas;
using namespace ambient;
//
// List of types T for which the p_dense_matrix<T> is tested
// (long long unsigned int causes problems in boost::iterator facade)
typedef boost::mpl::list<double> test_types;
typedef ambient::dim2 dim;

struct caveats {
     caveats(){ ambient::init();     }
    ~caveats(){ ambient::finalize(); }
};

BOOST_GLOBAL_FIXTURE( caveats );
/*
BOOST_AUTO_TEST_CASE_TEMPLATE( p_diag, T, test_types )
{
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1);
    p_diagonal_matrix<T> A(M_SIZE,1);

    ambient::push(ambient::one_l_scalapack_kernel, ambient::one_null_c_kernel, A.get_data() );
    ambient::playout();
    zout << A ;

    A.remove_rows(2,3);
    ambient::playout();
    zout << A ;

    A.resize(8,8);
    ambient::playout();
    A.resize(M_SIZE,M_SIZE);
    ambient::playout();
    zout << A ;
}
*/

BOOST_AUTO_TEST_CASE_TEMPLATE( p_diagonal_gemm, T, test_types )
{
    ambient::layout >> dim(2,2),dim(2,2),dim(2,2);
    p_diagonal_matrix<T> A(M_SIZE);
    p_dense_matrix<T>    B(M_SIZE,M_SIZE);
    p_dense_matrix<T>    C(M_SIZE,M_SIZE);

    blas::gemm(A,B,C);  
    ambient::playout();

    std::cout << A;
    std::cout << B;
    zout << "----\n";
    std::cout << C;
}


