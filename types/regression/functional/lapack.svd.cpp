
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

BOOST_AUTO_TEST_CASE_TEMPLATE( SVD_COMPARISON, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix pU(T::valuex,T::valuey);
    pMatrix pV(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sU(T::valuex,T::valuey);
    sMatrix sV(T::valuex,T::valuey);

    pDiagMatrix pS(T::valuex);
    sDiagMatrix sS((std::size_t)T::valuex); 

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    svd(pA,pU,pV,pS);
    svd(sA,sU,sV,sS);

    printf("Done %d x %d\n", (int)T::valuex, (int)T::valuey);
    BOOST_CHECK(sS == pS);
}
