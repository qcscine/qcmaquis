#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>
#include <cmath>

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"

#include "types/utils/matrix_cast.h"

#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( addition, T, test_types)
{
    ambient::layout >> dim(1,1), dim(1,1), dim(1,1);

    pMatrix pA(T::valuex,T::valuex);
    pMatrix pB(T::valuex,T::valuex);
    pDiagMatrix pC(T::valuex,T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sB(T::valuex,T::valuex);
    sDiagMatrix sC(T::valuex,T::valuex);

    maquis::types::algorithms::generate(sB,Rd); // Rd is rand generator static variable inside utilities
    sC.generate(Rd); // Rd is rand generator static variable inside utilities

    pB = maquis::traits::matrix_cast<pMatrix>(sB); // playout is inside the cast
    pC = maquis::traits::matrix_cast<pDiagMatrix>(sC); // playout is inside the cast
 
    gemm(pC,pB,pA);
    gemm(sC,sB,sA);

    

    BOOST_CHECK(pA==sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

