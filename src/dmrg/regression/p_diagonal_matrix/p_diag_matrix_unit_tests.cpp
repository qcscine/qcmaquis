#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_diagonal_matrix.h"

#include "p_dense_matrix/p_dense_matrix_algorithms.h"
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
/*
BOOST_AUTO_TEST_CASE_TEMPLATE( p_diag, T, test_types )
{
    ambient::layout >> dim3(1,1), dim3(2,2), dim3(10,1);
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

BOOST_AUTO_TEST_CASE_TEMPLATE( p_diag_gemm, T, test_types )
{
    ambient::layout >> dim3(1,1), dim3(2,2), dim3(10,1);
    p_dense_matrix<T> B(M_SIZE,M_SIZE);

    ambient::push(ambient::one_l_scalapack_kernel, ambient::one_null_c_kernel, B);
    ambient::playout();
    ambient::push(ambient::one_l_scalapack_kernel, ambient::one_null_c_kernel, B);
    ambient::playout();
   
    //zout << C  << std::endl;
  

//  blas::gemm(B,A,C);


  //  zout << A <<  std::endl;
  //  zout << B  << std::endl;
  //  std::cout << C;

}


