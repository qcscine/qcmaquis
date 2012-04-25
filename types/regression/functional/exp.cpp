
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/random.hpp>

#include "types/p_dense_matrix/p_dense_matrix.h"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/algorithms/algorithms.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"

#include "types/utils/matrix_cast.h"

#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( Exp, T, test_types)
{
    ambient::model >> dim(32,32), dim(32,32), dim(32,32);
    pDiagMatrix pA(T::valuex,0);
    sDiagMatrix sA(T::valuex,0);

    pA.get_data().set_init(ambient::random_i<typename T::dbl>);
    sA = maquis::traits::matrix_cast<sDiagMatrix>(pA); // playout is inside the cast

    sA = exp(sA);
    pA = exp(pA);

    BOOST_CHECK(pA==sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}
