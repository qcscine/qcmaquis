#include <cmath>
#include "utils/zout.hpp"
#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"
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

BOOST_AUTO_TEST_CASE_TEMPLATE( single_bench_gemm, T, test_types ) 
{
     ambient::layout >> dim(1,1), dim(1,1), dim(1,1);

     int argc = boost::unit_test::framework::master_test_suite().argc;
     char** argv = boost::unit_test::framework::master_test_suite().argv; 
     int NUM=std::atoi(argv[1]);
    
     ambient::p_dense_matrix<T> A(NUM,NUM);
     ambient::p_dense_matrix<T> B(NUM,NUM);
     ambient::p_dense_matrix<T> C_Ambient(NUM,NUM);
     ambient::p_dense_matrix<T> C_pBlas(NUM,NUM);

     A.set_init(ambient::random_i<T>);
     B.set_init(ambient::random_i<T>);
 
     Timer ta("Ambient_gemm.txt"); ta.begin();
     maquis::types::gemm(A,B,C_Ambient);
     ambient::playout();
     ta.end();

     Timer tb("Pblas_gemm.txt"); tb.begin();
//     maquis::types::pblas_gemm(A,B,C_pBlas);
//     ambient::playout();
     tb.end();     

     maquis::types::validation(C_Ambient,C_pBlas);
     ambient::playout();
     if(ambient::rank() == 0){
        save("single",ambient::size(),NUM, ta, tb );
     } 
}


BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_vector, T, test_types ) 
{
     ambient::layout >> dim(1,1), dim(1,1), dim(1,1); 
     char** argv = boost::unit_test::framework::master_test_suite().argv; 
     int argc = boost::unit_test::framework::master_test_suite().argc;
     int M, LENGTH;
     bool save_output;

     if(argc <= 1){ M = 4096; LENGTH = 8; save_output = false; }
     else{ M = std::atoi(argv[1]); LENGTH = std::atoi(argv[2]); save_output = true; }
     std::vector<ambient::p_dense_matrix<T> * > V;
     V.resize(LENGTH*4);

     for(int i = 0 ; i < LENGTH*4 ; i++){
         V[i] = (p_dense_matrix<T>*) new p_dense_matrix<T,MANUAL>(M,M);
         V[i]->set_init(ambient::random_i<T>);
     }

//     for(int i = 0 ; i < LENGTH ; i++) maquis::types::pblas_gemm(*V[i*4],*V[i*4+1],*V[i*4+3]);

     Timer ta("PBLAS series of GEMM"); ta.begin();
  //   ambient::playout();
     ta.end();

     for(int i = 0 ; i < LENGTH ; i++) maquis::types::gemm(*V[i*4],*V[i*4+1],*V[i*4+2]);

     Timer tb("Ambient series of GEMM"); tb.begin();
     ambient::playout();
     tb.end();

     for(int i = 0 ; i < LENGTH ; i++) maquis::types::validation(*V[i*4+2],*V[i*4+3]);
     ambient::playout();

     if(save_output && ambient::rank() == 0) save("vector.txt",ambient::size(),M, ta, tb );
}

