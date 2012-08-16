
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

BOOST_AUTO_TEST_CASE_TEMPLATE( addition, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix pQ(T::valuex,T::valuey);
    pMatrix pR(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sQ(T::valuex,T::valuey);
    sMatrix sR(T::valuex,T::valuey);


    pA.fill_random();
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
 
    qr(pA,pQ,pR);
    qr(sA,sQ,sR);

    BOOST_CHECK(sQ == pQ); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
    BOOST_CHECK(sR == pR); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}
