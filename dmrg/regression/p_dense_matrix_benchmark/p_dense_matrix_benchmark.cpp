#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "utils/timings.h"

#define BOOST_TEST_MODULE p_dense_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

using namespace blas;
using namespace ambient;

typedef boost::mpl::list<double> test_types;
typedef ambient::dim2 dim;

struct caveats {
    caveats(){ ambient::init();     }
   ~caveats(){ ambient::finalize(); }
};

BOOST_GLOBAL_FIXTURE( caveats );

BOOST_AUTO_TEST_CASE_TEMPLATE( single_gemm_test, T, test_types ) 
{ 
    ambient::layout >> dim(1,1), dim(1,1), dim(1,1); 
    int M = 2;
    int N = 2;

    p_dense_matrix<T> V1(M,N*2); 
    p_dense_matrix<T> V2(N*2,N*2); 
    p_dense_matrix<T> V3(M,N*2); 
    p_dense_matrix<T> V4(M,N*2); 
    blas::gemm(V1,V2,V3);
    blas::pblas_gemm(V1,V2,V4);
    blas::validation(V3,V4);

    Timer b("Single GEMM: Ambient; PDGEMM;");
    b.begin(); 
    ambient::playout();
    b.end();

    //if(ambient::rank() == 0) b.save(ambient::size(),M);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_vector, T, test_types ) 
{
     ambient::layout >> dim(1,1), dim(1,1), dim(1,1); 
     int LENGTH = 8;
     int M = 128;

     std::vector<ambient::p_dense_matrix<T> * > V;
     V.resize(LENGTH*4);

     for(int i = 0 ; i < LENGTH*4 ; i++) V[i] = (p_dense_matrix<T>*) new p_dense_matrix<T,MANUAL>(M,M);

     for(int i = 0 ; i < LENGTH ; i++) blas::pblas_gemm(*V[i*4],*V[i*4+1],*V[i*4+3]);
     Timer tp("PBLAS series of GEMM"); tp.begin();
     ambient::playout();
     tp.end();
     //if(ambient::rank() == 0) ta.save(ambient::size(),M);
     for(int i = 0 ; i < LENGTH ; i++) blas::gemm(*V[i*4],*V[i*4+1],*V[i*4+2]);
     Timer ta("Ambient series of GEMM"); ta.begin();
     ambient::playout();
     ta.end();

     for(int i = 0 ; i < LENGTH ; i++) blas::validation(*V[i*4+2],*V[i*4+3]);
     ambient::playout();
}
