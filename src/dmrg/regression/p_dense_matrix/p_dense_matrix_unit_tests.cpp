#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.hpp"
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
typedef ambient::dim2 dim;

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
    ambient::layout >> dim(1,1), dim(1,1), dim(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);
    p_dense_matrix<T> d(M_SIZE,M_SIZE);

    c = a + b;
    try{
    c(5,5) = 13.0;
    c(6,5) = 14.0;
    }catch(...){}
    a.remove_rows(128,128);
    b.remove_rows(128,128);
    c.remove_rows(128,128);
    c = a + b;
    a.resize(640,512);
    b.resize(640,512);
    c.resize(640,512);
    c = a + b;
    ambient::playout();
}


BOOST_AUTO_TEST_CASE_TEMPLATE( sql_test, T, test_types )
{
    ambient::layout >> dim(1,1), dim(1,1), dim(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);

    c = a + b;
    c = a + b;
    ambient::playout();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( print_test, T, test_types )
{
    ambient::layout >> dim(1,1), dim(1,1), dim(10,1);

    p_dense_matrix<T> a(8,8);
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,a); 
    a.remove_rows(0,1);

    std::cout << a;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SVD_test, T, test_types ) 
{ 
    ambient::layout >> dim(1,1), dim(1,1), dim(10,1); 
 
    p_dense_matrix<T> A(M_SIZE,M_SIZE); 
    p_dense_matrix<T> U(M_SIZE,M_SIZE); 
    p_dense_matrix<T> V(M_SIZE,M_SIZE); 
 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,A); 
 
    typename::associated_diagonal_matrix<p_dense_matrix<T> >::type S; 
 
 
    blas::svd(A,U,V,S); 
 
    ambient::playout(); 
//    std::cout << S ; 
 
}*/


/*BOOST_AUTO_TEST_CASE_TEMPLATE( heap_manual_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    p_dense_matrix<T,ambient::MANUAL>* A = new p_dense_matrix<T, ambient::MANUAL>(M_SIZE,M_SIZE);
    p_dense_matrix<T,ambient::MANUAL>* U = new p_dense_matrix<T, ambient::MANUAL>(M_SIZE,M_SIZE); 
    p_dense_matrix<T,ambient::MANUAL>* V = new p_dense_matrix<T, ambient::MANUAL>(M_SIZE,M_SIZE); 
 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,*A); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,*U); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,*V); 
    *A = *U + *V;
    ambient::playout();
    delete A; delete U; delete V;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( heap_auto_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    p_dense_matrix<T,ambient::WEAK>* A = new p_dense_matrix<T,ambient::WEAK>(M_SIZE,M_SIZE);
    p_dense_matrix<T,ambient::WEAK>* U = new p_dense_matrix<T,ambient::WEAK>(M_SIZE,M_SIZE); 
    p_dense_matrix<T,ambient::WEAK>* V = new p_dense_matrix<T,ambient::WEAK>(M_SIZE,M_SIZE); 

    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,*A); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,*U); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,*V); 
    *A = *U + *V;
    ambient::playout();
    delete A;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( stack_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    p_dense_matrix<T> A(M_SIZE,M_SIZE);
    p_dense_matrix<T> U(M_SIZE,M_SIZE); 
    p_dense_matrix<T> V(M_SIZE,M_SIZE); 

    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,A); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,U); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,V); 
    A = U + V;
    ambient::playout();
}*/

BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 

    size_t task_size = M_SIZE;

    p_dense_matrix<T> A(task_size,task_size);
    p_dense_matrix<T> B(task_size,task_size); 
    p_dense_matrix<T> C(task_size,task_size);

    p_dense_matrix<T> A2(task_size,task_size);
    p_dense_matrix<T> B2(task_size,task_size); 
    p_dense_matrix<T> C2(task_size,task_size); 

    C = A * B;
    ambient::push(ambient::check_gemm_l_kernel,ambient::check_gemm_c_kernel,C); 
    C2 = A2 * B2;
    ambient::push(ambient::check_gemm_l_kernel,ambient::check_gemm_c_kernel,C2); 
    ambient::playout();
}
