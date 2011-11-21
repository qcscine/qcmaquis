#include <cmath>
#include "utils/zout.hpp"
#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "utils/timings.h"

#define T double

using namespace blas;
using namespace ambient;

typedef ambient::dim2 dim;

int main(int argc, char* argv[])
{
     ambient::init();
     ambient::layout >> dim(1,1), dim(1,1), dim(1,1); 
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

//     for(int i = 0 ; i < LENGTH ; i++)
//           maquis::types::pblas_gemm(*V[i*4],*V[i*4+1],*V[i*4+3]);

     Timer ta("PBLAS series of GEMM");
     ta.begin();
  //   ambient::playout();
     ta.end();

     for(int i = 0 ; i < LENGTH ; i++)
         maquis::types::gemm(*V[i*4],*V[i*4+1],*V[i*4+2]);

     Timer tb("Ambient series of GEMM");
     tb.begin();
     ambient::playout();
     tb.end();

//     for(int i = 0 ; i < LENGTH ; i++) 
//         maquis::types::validation(*V[i*4+2],*V[i*4+3]);
//     ambient::playout();

     if(save_output && ambient::rank() == 0) save("vector.txt",ambient::size(),M, ta, tb );

     ambient::finalize();
}

