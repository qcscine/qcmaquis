#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "p_dense_matrix/concept/matrix_interface.hpp"
#include "p_dense_matrix/concept/resizable_matrix_interface.hpp"
#include "utils/timings.h"

#define BOOST_TEST_MODULE p_dense_matrix
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <complex>
#include <numeric>

#define M_SIZE 128
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


/*
BOOST_AUTO_TEST_CASE_TEMPLATE(VALIDATION, T, test_types ) 
{ 
    ambient::layout >> dim3(10,5), dim3(2,2), dim3(10,1); 

//    int M = atoi(boost::unit_test::framework::master_test_suite().argv[1]); 
     int  M = 256 ;
    p_dense_matrix<T> A(M,M); 

//    p_dense_matrix<T> B(M,M); 
//    p_dense_matrix<T> C(M,M); 
//    p_dense_matrix<T> D(M,M); 
 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,A); 
//    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,B); 
//    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,C); 
//   ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,D); 
    
//    blas::gemm(A,B,C);
//    blas::ambient_gemm(A,B,D);
  //  ambient::push(ambient::validation_l_kernel, ambient::validation_c_kernel, D,C );
    ambient::playout(); 
    //std::cout << C << std::endl;
    //std::cout << D << std::endl;
}
*/

BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_AMBIENT_test2, T, test_types ) 
{ 

    ambient::layout >> dim3(10,5), dim3(8,8), dim3(10,1); 
   
    //int M = atoi(boost::unit_test::framework::master_test_suite().argv[1]); 
    int M = 1024;

    zout << M ;  
  /* 
    p_dense_matrix<T> A1(M,M); 
    p_dense_matrix<T> A2(M,M); 
    p_dense_matrix<T> A3(M,M); 

    M  *= 2;
    p_dense_matrix<T> B1(M,M); 
    p_dense_matrix<T> B2(M,M); 
    p_dense_matrix<T> B3(M,M); 

    M  *= 2;
    p_dense_matrix<T> C1(M,M); 
    p_dense_matrix<T> C2(M,M); 
    p_dense_matrix<T> C3(M,M); 
*/
    M  *= 2;
    p_dense_matrix<T> D1(M,M); 
    p_dense_matrix<T> D2(M,M); 
    p_dense_matrix<T> D3(M,M); 
/*
    M  /= 8;
    M *= 6;
    p_dense_matrix<T> E1(M,M); 
    p_dense_matrix<T> E2(M,M); 
    p_dense_matrix<T> E3(M,M); 
 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,A1); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,A2); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,A3); 
  
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,B1); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,B2); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,B3); 

    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,C1); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,C2); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,C3); 

    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,D1); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,D2); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,D3); 

    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,E1); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,E2); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,E3); 
*/
//    blas::gemm(A1,A2,A3);
//    blas::gemm(B1,B2,B3);
//    blas::gemm(C1,C2,C3);
    blas::gemm(D1,D2,D3);
/*
    blas::gemm(E1,E2,E3);
*/

 //   blas::ambient_gemm(C1,c2,C3);
 //   blas::ambient_gemm(D1,D2,D3);
 //   blas::ambient_gemm(E1,E2,E3);

    Timer b("AMBIENT_main");
    b.begin(); 
    ambient::playout();
    b.end();
 
    if(ambient::rank() == 0)
    {
	b.save(ambient::size(),M);
    }   
  /*
   ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,AV3); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,BV3); 

    blas::gemm(A1,A2,AV3);
    blas::gemm(B1,B2,BV3);

    ambient::push(ambient::validation_l_kernel, ambient::validation_c_kernel, A3,AV3 );
    ambient::push(ambient::validation_l_kernel, ambient::validation_c_kernel, B3,BV3 );

    ambient::playout();
//    if(ambient::rank() == 0)
    {
	b.save(ambient::size(),M);
    }   
*/
}

