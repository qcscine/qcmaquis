#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>
#include <cmath>
#include <limits>

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "types/p_dense_matrix/algorithms.hpp"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"

#include "types/utils/matrix_cast.h"

#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( addition, T, test_types)
{
    ambient::layout >> dim(1,1), dim(1,1), dim(1,1);

    pMatrix pA(T::valuex,T::valuex);
    sMatrix sA(T::valuex,T::valuex);

    typename T::dbl a,b;

    pA.set_init(ambient::random_i<typename T::dbl>);
    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast

    a = trace(pA); 
    b = trace(sA); 
 
    ambient::playout();
    BOOST_CHECK_CLOSE( a, b,0.000000000001); 
}


