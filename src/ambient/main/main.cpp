#include <mpi.h>
#include <complex>
#include <iostream>
#include <cmath>
#include "utils/zout.hpp"

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"

#include "types/utils/matrix_cast.h" 

using namespace blas;
using namespace ambient;

typedef ambient::dim2 dim;

//#define  T double
#define  T std::complex<double

int main(int argc, char* argv[])
{
    ambient::init();

     ambient::layout >> dim(1,1), dim(1,1), dim(1,1);

     int NUM=std::atoi(argv[1]);
    
     p_dense_matrix<T> pA(NUM,NUM);
     pA.set_init(ambient::random_i<T>);
     std::cout << pA << std::endl;
  //   ambient::playout();   
  //   ambient::playout();   
/*
    // pA += pC;

         std::cout << "0 OK " << std::endl;

         std::cout << "1  OK " << std::endl;
     B = matrix_cast<ambient::dense_matrix<T> >(pB);

         std::cout << "2  OK " << std::endl;
     C = matrix_cast<ambient::dense_matrix<T> >(pC);

     maquis::types::gemm(A,B,D);

     if(C == D)
         std::cout << " OK " << std::endl;
     ambient::playout();  
     zout << C << std::endl;
*/
     ambient::finalize();
}

