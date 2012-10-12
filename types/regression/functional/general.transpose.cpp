#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( TRANSPOSE, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    transpose_inplace(pA); 
    transpose_inplace(sA); 

    BOOST_CHECK(pA==sA);
}

