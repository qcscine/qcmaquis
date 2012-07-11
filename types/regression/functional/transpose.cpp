#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( Transpose, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix ptA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);
    sMatrix stA(T::valuex,T::valuey);

    pA.fill_random();

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast

    ptA = transpose(pA); 
    stA = transpose(sA); 

    BOOST_CHECK(ptA==stA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

