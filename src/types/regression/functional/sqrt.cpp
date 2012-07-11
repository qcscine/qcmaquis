#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/random.hpp>

#include "alps/numeric/matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( Sqrt, T, test_types)
{
    pDiagMatrix pA(T::valuex,0);
    sDiagMatrix sA(T::valuex,0);

    pA.get_data().fill_random();
    sA = maquis::traits::matrix_cast<sDiagMatrix>(pA); // playout is inside the cast

    sA = sqrt(sA);
    pA = sqrt(pA);

    BOOST_CHECK(pA==sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}
