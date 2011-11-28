#include <mpi.h>

#include <iostream>
#include <cmath>
//#include "utils/zout.hpp"
#include "ambient/ambient.h"
#include "utils/zout.hpp"

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/concept/matrix_interface.hpp"
#include "types/p_dense_matrix/algorithms.hpp"
#include "types/p_dense_matrix/concept/resizable_matrix_interface.hpp"
#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "types/p_dense_matrix/p_diagonal_matrix.h"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/dense_matrix_algorithms.h"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#include "types/utils/matrix_cast.h"

#define  T double

//using namespace maquis::types::algorithms;
using namespace maquis::types;
using namespace ambient;

typedef maquis::types::dense_matrix<T> sMatrix;
typedef maquis::types::diagonal_matrix<T> sdMatrix;
typedef maquis::types::p_dense_matrix<T> pMatrix;
typedef maquis::types::p_diagonal_matrix<T> pdMatrix;

typedef ambient::dim2 dim;

namespace Random{
   struct random {
       random(){};
       double operator()(){return drand48();}
       int IntRd(){return mrand48();}
   };
}

static Random::random Rd; //one generator is enough so static

int main(int argc, char* argv[])
{
    ambient::init();

    ambient::layout >> dim(1,1), dim(1,1), dim(1,1);

    int NUM=std::atoi(argv[1]);
    
    pMatrix pA(NUM,NUM);
    pMatrix pB(NUM,NUM);
    pMatrix pC(NUM,NUM);

    sMatrix sA(NUM,NUM);
    sMatrix sB(NUM,NUM);
    sMatrix sC(NUM,NUM);

    pA.set_init(ambient::random_i<double>);
    pB.set_init(ambient::random_i<double>);
    pC.set_init(ambient::null_i<double>);

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
    sB = maquis::traits::matrix_cast<sMatrix>(pB); // playout is inside the cast
//    sC = maquis::traits::matrix_cast<sMatrix>(pC); // playout is inside the cast
 
    sC = sA + sB; 
    pC = pA + pB; 
//    zout << sA << std::endl;  
    zout << " --------sA------------ " << std::endl;
    zout << sA << std::endl;  
    zout << " --------sB------------ " << std::endl;
    zout << sB << std::endl;  
    zout << " --------sC------------ " << std::endl;
    zout << sC << std::endl;  
    zout << " ---------------------- " << std::endl;
    zout << " --------pA------------ " << std::endl;
    std::cout << pA << std::endl;  
    zout << " --------pB------------ " << std::endl;
    std::cout << pB << std::endl;  
    zout << " --------pC------------ " << std::endl;
    std::cout << pC << std::endl;  
//    zout << " ---------------------- " << std::endl;
//    std::cout << pC << std::endl;  
/*   
    if( pC == sC )
        zout << " ok "   
 */
    ambient::finalize();

}
