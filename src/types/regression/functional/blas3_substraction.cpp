#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( substraction, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix pB(T::valuex,T::valuey);
    pMatrix pC(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sB(T::valuex,T::valuey);
    sMatrix sC(T::valuex,T::valuey);

    generate(pA);
    generate(pB);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    sB = maquis::bindings::matrix_cast<sMatrix>(pB);
    sC = maquis::bindings::matrix_cast<sMatrix>(pC);
 
    sC = sA - sB; 
    pC = pA - pB; 

    BOOST_CHECK(pC == sC); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

BOOST_AUTO_TEST_CASE_TEMPLATE( substraction_assign, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix pB(T::valuex,T::valuey);
    pMatrix pC(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sB(T::valuex,T::valuey);
    sMatrix sC(T::valuex,T::valuey);

    generate(pA);
    generate(pB);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    sB = maquis::bindings::matrix_cast<sMatrix>(pB);
    sC = maquis::bindings::matrix_cast<sMatrix>(pC);
 
    sA -= sB; 
    pA -= pB; 

    BOOST_CHECK(pA == sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

