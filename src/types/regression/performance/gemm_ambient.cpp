#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>
#include <cmath>

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/algorithms.hpp"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/utils/matrix_cast.h"

#include "utils/timings.h"
#include "utilities.h"

typedef ambient::dim2 dim;
typedef maquis::types::dense_matrix<ValueType> sMatrix;
typedef maquis::types::p_dense_matrix<ValueType> pMatrix;
typedef maquis::types::diagonal_matrix<ValueType> sDiagMatrix;
typedef maquis::types::p_diagonal_matrix<ValueType> pDiagMatrix;

BOOST_AUTO_TEST_CASE_TEMPLATE( test, T, test_types)
{
    ambient::model >> dim(256,256), dim(256,256), dim(256,256);

    pMatrix pA(T::ValueX,T::ValueY);
    pMatrix pB(T::ValueX,T::ValueY);
    pMatrix pC(T::ValueX,T::ValueY);

    sMatrix sA(T::ValueX,T::ValueY);
    sMatrix sB(T::ValueX,T::ValueY);

    pA.set_init(ambient::random_i<typename T::value_type>);
    pB.set_init(ambient::random_i<typename T::value_type>);

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
    sB = maquis::traits::matrix_cast<sMatrix>(pB); // playout is inside the cast
    ambient::playout();

    maquis::types::algorithms::gemm(pA,pB,pC); 

    TimerPTH t1(" ambient");
    t1.begin();
    ambient::playout();
    t1.end();

    double time = t1.GetTime();
    double gfl  = GFlopsGemm(T::ValueX,T::ValueX,T::ValueX,time);
    save("TimeGemmAmbient.txt",t1,gfl,T::ValueX,T::ValueY,T::ValueThread); 
}

