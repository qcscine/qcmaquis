#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "types/p_dense_matrix/p_dense_matrix.h"

#include "alps/numeric/matrix.hpp"

#include "types/utils/bindings.hpp"

#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( resize, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    pA.fill_random();
    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast

    sA.resize(4,2);
    pA.resize(4,2);

    BOOST_CHECK(pA == sA);
}


