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

#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( dgemm, T, test_types)
{
    ambient::model >> dim(128,128), dim(128,128), dim(128,128);

    pMatrix pA(T::valuex,T::valuex);
    pMatrix pB(T::valuex,T::valuex);
    pMatrix pC(T::valuex,T::valuex);


    sMatrix sA(T::valuex,T::valuex);
    sMatrix sB(T::valuex,T::valuex);
    sMatrix sC(T::valuex,T::valuex);

    pA.set_init(ambient::random_i<typename T::dbl>);
    pB.set_init(ambient::random_i<typename T::dbl>);

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
    sB = maquis::traits::matrix_cast<sMatrix>(pB); // playout is inside the cast

    printf("BEGIN OF GEMM!\n"); 
    maquis::types::algorithms::gemm(pA,pB,pC);
    ambient::playout();
    printf("END OF GEMM!\n"); 
    printf("BEGIN OF SERIAL GEMM!\n"); 
    maquis::types::gemm(sA,sB,sC); // to fix 
    printf("END OF SERIAL GEMM!\n"); 
    BOOST_CHECK(pC==sC); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

