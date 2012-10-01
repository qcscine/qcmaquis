#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( dgemm, T, test_types)
{
    pMatrix pA(T::valuex,T::valuex);
    pMatrix pB(T::valuex,T::valuex);
    pMatrix pC(T::valuex,T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sB(T::valuex,T::valuex);
    sMatrix sC(T::valuex,T::valuex);

    generate(pA);
    generate(pB);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    sB = maquis::bindings::matrix_cast<sMatrix>(pB);

    ambient::numeric::gemm(pA,pB,pC);
    ambient::sync();
    gemm(sA,sB,sC);
    BOOST_CHECK(pC==sC);
}

