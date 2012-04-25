
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/algorithms.hpp"
#include "types/p_dense_matrix/p_diagonal_matrix.h"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/dense_matrix_algorithms.h"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#include "types/utils/matrix_cast.h"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( syev, T, test_types)
{
    ambient::model >> dim(32,32), dim(32,32), dim(32,32);

    pMatrix pA(T::valuex,T::valuex);
    pMatrix pV(T::valuex,T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sV(T::valuex,T::valuex);

    pA.set_init(ambient::random_i<typename T::dbl>);
    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
 
    typename maquis::types::associated_diagonal_matrix<pMatrix>::type pE; 
    typename maquis::types::associated_diagonal_matrix<sMatrix>::type sE;

    pE.resize(T::valuex,T::valuex); 
    sE.resize(T::valuex,T::valuex); 
 
    maquis::types::syev(pA,pV,pE); // to modify the algo we need the reverse inside !
    maquis::types::algorithms::syev(sA,sV,sE);
     
    maquis::types::diagonal_matrix<typename T::dbl> sE2(maquis::traits::matrix_cast<maquis::types::diagonal_matrix<typename T::dbl> >(pE));
    BOOST_CHECK(sE == pE);
}
