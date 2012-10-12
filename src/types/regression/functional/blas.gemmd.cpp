#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( GEMM_DIAGONAL, T, test_types)
{
    pMatrix pA(T::valuex,T::valuex);
    pMatrix pB(T::valuex,T::valuex);
    pDiagMatrix pC(T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sB(T::valuex,T::valuex);
    sDiagMatrix sC((std::size_t)T::valuex);

    generate(pB);
    generate(pC);

    sB = maquis::bindings::matrix_cast<sMatrix>(pB);
    sC = maquis::bindings::matrix_cast<sDiagMatrix>(pC);

    ambient::numeric::gemm(pC,pB,pA);
    gemm(sC,sB,sA);

    BOOST_CHECK(pA==sA);
}

