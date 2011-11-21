
#include <cmath>
#include "utils/zout.hpp"
#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "utils/timings.h"

using namespace blas;
using namespace ambient;

typedef ambient::dim2 dim;

#define  T double

int main(int argc, char* argv[])
{
    ambient::init();

     ambient::layout >> dim(1,1), dim(1,1), dim(1,1);

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
 //    ambient::playout();
     tb.end();     

//     maquis::types::validation(C_Ambient,C_pBlas);
//     ambient::playout();
     if(ambient::rank() == 0){
        save("single",ambient::size(),NUM, ta, tb );
     } 

     ambient::finalize();
}

