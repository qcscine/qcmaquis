
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "ambient/numeric/matrix.hpp"
#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( LQ_COMPARISON, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix pL(T::valuex,T::valuey);
    pMatrix pQ(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sL(T::valuex,T::valuey);
    sMatrix sQ(T::valuex,T::valuey);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
 
    lq(pA,pL,pQ);
    lq(sA,sL,sQ);

    BOOST_CHECK(sL == pL);
    BOOST_CHECK(sQ == pQ);
}
