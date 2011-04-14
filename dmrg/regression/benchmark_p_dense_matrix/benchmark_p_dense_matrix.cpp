
#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.h"
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


/*BOOST_AUTO_TEST_CASE_TEMPLATE( summ_operation_test, T, test_types )
{
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1);

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
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);

    c = a + b;
    c = a + b;
    ambient::playout();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( print_test, T, test_types )
{
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1);

    p_dense_matrix<T> a(8,8);
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,a); 
    a.remove_rows(0,1);

    std::cout << a;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SVD_test, T, test_types ) 
{ 
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1); 
 
    p_dense_matrix<T> A(M_SIZE,M_SIZE); 
    p_dense_matrix<T> U(M_SIZE,M_SIZE); 
    p_dense_matrix<T> V(M_SIZE,M_SIZE); 
 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,A); 
 
    typename::associated_diagonal_matrix<p_dense_matrix<T> >::type S; 
 
    blas::svd(A,U,V,S); 

    ambient::playout();
}


BOOST_AUTO_TEST_CASE_TEMPLATE(VALIDATION, T, test_types ) 
{ 
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1); 

    int M = atoi(boost::unit_test::framework::master_test_suite().argv[1]); 
    M = 256 ;
    p_dense_matrix<T> A(M,M); 
    p_dense_matrix<T> B(M,M); 
    p_dense_matrix<T> C(M,M); 
    p_dense_matrix<T> D(M,M); 
 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,A); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,B); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,C); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,D); 
    
    blas::gemm(A,B,C);
//    blas::ambient_gemm(A,B,D);
  //  ambient::push(ambient::validation_l_kernel, ambient::validation_c_kernel, D,C );
    ambient::playout(); 
    //std::cout << C << std::endl;
    //std::cout << D << std::endl;
}


<<<<<<< .mine
BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_AMBIENT, T, test_types ) 
=======
BOOST_AUTO_TEST_CASE_TEMPLATE( heap_manual_test, T, test_types ) 
>>>>>>> .r637
{ 

    Timer b("AMBIENT");
    b.begin(); 
    ambient::layout >> dim3(10,5), dim3(1,1), dim3(10,1); 

    int M = atoi(boost::unit_test::framework::master_test_suite().argv[1]); 
    p_dense_matrix<T> A(M,M); 
    p_dense_matrix<T> U(M,M); 
    p_dense_matrix<T> V(M,M); 
<<<<<<< .mine
 
=======
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
    ambient::layout >> dim3(10,5), dim3(2,2), dim3(10,1); 
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
    ambient::layout >> dim3(10,5), dim3(2,2), dim3(10,1); 
    p_dense_matrix<T> A(M_SIZE,M_SIZE);
    p_dense_matrix<T> U(M_SIZE,M_SIZE); 
    p_dense_matrix<T> V(M_SIZE,M_SIZE); 

>>>>>>> .r637
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,A); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,U); 
<<<<<<< .mine
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,V); 
    
    blas::ambient_gemm(A,U,V);
=======
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,V); 
    A = U + V;
>>>>>>> .r637
    ambient::playout();
    b.end();
     
    if(ambient::rank() == 0)
    {
	b.save(ambient::size(),M);
    } 
  
}
*/

BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_AMBIENT_test2, T, test_types ) 
{ 

    Timer b("AMBIENT");
    b.begin(); 
    ambient::layout >> dim3(10,5), dim3(2,2), dim3(10,1); 
   
    //int M = atoi(boost::unit_test::framework::master_test_suite().argv[1]); 
    int M = 256;

    zout << M ;  
    
    p_dense_matrix<T> A1(M,M); 
    p_dense_matrix<T> A2(M,M); 
    p_dense_matrix<T> A3(M,M); 
    p_dense_matrix<T> AV3(M,M); 

    M  *= 2;
    p_dense_matrix<T> B1(M,M); 
    p_dense_matrix<T> B2(M,M); 
    p_dense_matrix<T> B3(M,M); 
    p_dense_matrix<T> BV3(M,M); 

    M  *= 2;
    p_dense_matrix<T> C1(M,M); 
    p_dense_matrix<T> C2(M,M); 
    p_dense_matrix<T> C3(M,M); 

    M  *= 2;
    p_dense_matrix<T> D1(M,M); 
    p_dense_matrix<T> D2(M,M); 
    p_dense_matrix<T> D3(M,M); 
/*
    M  *= 2;
    p_dense_matrix<T> E1(M,M); 
    p_dense_matrix<T> E2(M,M); 
    p_dense_matrix<T> E3(M,M); 
*/
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
/*
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,E1); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel,E2); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,E3); 
*/
    blas::ambient_gemm(A1,A2,A3);
    blas::ambient_gemm(B1,B2,B3);
    blas::ambient_gemm(C1,C2,C3);
    blas::ambient_gemm(D1,D2,D3);
 //   blas::ambient_gemm(E1,E2,E3);


 //   blas::ambient_gemm(C1,c2,C3);
 //   blas::ambient_gemm(D1,D2,D3);
 //   blas::ambient_gemm(E1,E2,E3);

    ambient::playout();
    b.end();
/* 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,AV3); 
    ambient::push(ambient::init_double_l_kernel,ambient::init_double_c_kernel_0,BV3); 

    blas::gemm(A1,A2,AV3);
    blas::gemm(B1,B2,BV3);

    ambient::push(ambient::validation_l_kernel, ambient::validation_c_kernel, A3,AV3 );
    ambient::push(ambient::validation_l_kernel, ambient::validation_c_kernel, B3,BV3 );

    ambient::playout();
*/
//    if(ambient::rank() == 0)
    {
	b.save(ambient::size(),M);
    }   
}



