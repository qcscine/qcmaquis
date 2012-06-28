#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "types/utils/bindings.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( dgemm, T, test_types)
{
    pMatrix pA(T::valuex,T::valuex);
    pMatrix pB(T::valuex,T::valuex);
    pMatrix pC(T::valuex,T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sB(T::valuex,T::valuex);
    sMatrix sC(T::valuex,T::valuex);

    pA.fill_random();
    pB.fill_random();

    sA = maquis::traits::matrix_cast<sMatrix>(pA);
    sB = maquis::traits::matrix_cast<sMatrix>(pB);

    using maquis::types::NoTranspose;
    maquis::types::gemm<NoTranspose,NoTranspose>(pA,pB,pC);
    ambient::playout();
    gemm(sA,sB,sC);
    BOOST_CHECK(pC==sC); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

