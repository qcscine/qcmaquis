
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "types/p_dense_matrix/p_dense_matrix.h"

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"

#include "types/utils/bindings.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( syev_comparison, T, test_types)
{
    pMatrix pA(T::valuex,T::valuex);
    pMatrix pV(T::valuex,T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sV(T::valuex,T::valuex);

    pA.fill_random();
    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
 
    typename alps::numeric::associated_diagonal_matrix<pMatrix>::type pE(T::valuex,T::valuex); 
    typename alps::numeric::associated_diagonal_matrix<sMatrix>::type sE(T::valuex,T::valuex);

    syev(pA,pV,pE); // to modify the algo we need the reverse inside !
    syev(sA,sV,sE);
     
    BOOST_CHECK(sE == pE);
}
