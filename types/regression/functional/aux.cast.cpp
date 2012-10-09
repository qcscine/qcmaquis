
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/random.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( cast_p2s, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);
    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( cast_s2p, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);
    generate(sA,Rd); // Rd is rand generator static variable inside utilities
    pA = maquis::bindings::matrix_cast<pMatrix>(sA);
    BOOST_CHECK(sA==pA);
}
/*
BOOST_AUTO_TEST_CASE_TEMPLATE( cast_p2s_diag, T, test_types)
{
    pDiagMatrix pA(T::valuex);
    sDiagMatrix sA((std::size_t)T::valuex);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sDiagMatrix>(pA);
    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( cast_s2p_diag, T, test_types)
{
    pDiagMatrix pA(T::valuex);
    sDiagMatrix sA((std::size_t)T::valuex);
   
    sA.generate(Rd); // Rd is rand generator static variable inside utilities
    pA = maquis::bindings::matrix_cast<pDiagMatrix>(sA);
    BOOST_CHECK(sA==pA);
}
*/
