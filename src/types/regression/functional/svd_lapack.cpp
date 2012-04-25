
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "types/p_dense_matrix/p_dense_matrix.h"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/dense_matrix_algorithms.h"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/algorithms/algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#include "types/utils/matrix_cast.h"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( addition, T, test_types)
{
    ambient::model >> dim(32,32), dim(32,32), dim(32,32);

    pMatrix pA(T::valuex,T::valuey);
    pMatrix pU(T::valuex,T::valuey);
    pMatrix pV(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sU(T::valuex,T::valuey);
    sMatrix sV(T::valuex,T::valuey);

    pA.set_init(ambient::random_i<typename T::dbl>);
    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast

    typename maquis::types::associated_diagonal_matrix<pMatrix>::type pS;
    typename maquis::types::associated_diagonal_matrix<sMatrix>::type sS; 
 
    maquis::types::svd(pA,pU,pV,pS);
    maquis::types::algorithms::svd(sA,sU,sV,sS);
  
    BOOST_CHECK(sS == pS); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}
