#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "types/p_dense_matrix/p_dense_matrix.h"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"

#include "types/utils/bindings.hpp"

#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( addition, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix pB(T::valuex,T::valuey);
    pMatrix pC(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sB(T::valuex,T::valuey);
    sMatrix sC(T::valuex,T::valuey);

    pA.fill_random();
    pB.fill_random();

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
    sB = maquis::traits::matrix_cast<sMatrix>(pB); // playout is inside the cast
    sC = maquis::traits::matrix_cast<sMatrix>(pC); // playout is inside the cast

    sC = sA + sB;
    pC = pA + pB; 

    BOOST_CHECK(pC == sC); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

BOOST_AUTO_TEST_CASE_TEMPLATE( addition_assign, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix pB(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sB(T::valuex,T::valuey);

    pA.fill_random();
    pB.fill_random();

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
    sB = maquis::traits::matrix_cast<sMatrix>(pB); // playout is inside the cast
 
    sA += sB; 
    pA += pB; 

    BOOST_CHECK(pA == sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

