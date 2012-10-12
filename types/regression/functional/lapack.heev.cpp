
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( HEEV_COMPARISON, T, test_types)
{
    pMatrix pA(T::valuex,T::valuex);
    pMatrix pV(T::valuex,T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sV(T::valuex,T::valuex);

    pDiagMatrix pE(T::valuex); 
    sDiagMatrix sE((std::size_t)T::valuex);
 
    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    heev(pA,pV,pE);
    heev(sA,sV,sE);
     
    BOOST_CHECK(pE == sE);
}
