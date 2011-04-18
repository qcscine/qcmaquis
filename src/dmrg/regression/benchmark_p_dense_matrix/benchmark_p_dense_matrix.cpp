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
#include <vector>

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


//G
// List of type pairs <T,U> for which the mixed type matrix vector multiplication is tested.
//
typedef boost::mpl::list<type_pairs::IntDouble, type_pairs::DoubleInt> test_type_pairs;

BOOST_GLOBAL_FIXTURE( AmbientConfig );


BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_ambient, T, test_types ) 
{ 
//    ambient::init();

    ambient::layout >> dim3(10,5), dim3(8,8), dim3(10,1); 
   
    int M = 1024;

    M  *= 2;
    p_dense_matrix<T> D1(M,M); 
    p_dense_matrix<T> D2(M,M); 
    p_dense_matrix<T> D3(M,M); 
    p_dense_matrix<T> Ds3(M,M); 

    blas::gemm(D1,D2,D3);
    blas::pblas_gemm(D1,D2,Ds3);
    blas::validation(D3,Ds3);
   
    Timer b("AMBIENT_main");
    b.begin(); 
    ambient::playout();
    b.end();
 
    if(ambient::rank() == 0)
    {
	b.save(ambient::size(),M);
    }   

//    ambient::finalize();

}


BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_vector, T, test_types ) 
{
//   ambient::init();

     ambient::layout >> dim3(10,5), dim3(8,8), dim3(10,1); 
     int SIZE = 4;
     int M = 2048;

     std::vector<ambient::p_dense_matrix<T> * > V1;
     std::vector<ambient::p_dense_matrix<T> * > V2;
     std::vector<ambient::p_dense_matrix<T> * > V3;
     std::vector<ambient::p_dense_matrix<T> * > V4;

     V1.resize(SIZE);
     V2.resize(SIZE);
     V3.resize(SIZE);
     V4.resize(SIZE);

     for(int i = 0 ; i < SIZE ; i++)
     {
        V1[i] = new ambient::p_dense_matrix<T>(M,M);
        V2[i] = new ambient::p_dense_matrix<T>(M,M);
        V3[i] = new ambient::p_dense_matrix<T>(M,M);
        V4[i] = new ambient::p_dense_matrix<T>(M,M);
     }

     for(int i = 0 ; i < SIZE ; i++)
     {
        blas::gemm(*V1[i],*V2[i],*V3[i]);
     }
     
     Timer t("Test_vector_gemm");
     t.begin();
     ambient::playout();
     t.end();

     if(ambient::rank() == 0)
     {
 	 t.save(ambient::size(),M);
     }   

     for(int i = 0 ; i < SIZE ; i++)
     {
        blas::pblas_gemm(*V1[i],*V2[i],*V4[i]);
        blas::validation(*V3[i],*V4[i]);
     }

     ambient::playout();
 
//     ambient::finalize(); 

}



