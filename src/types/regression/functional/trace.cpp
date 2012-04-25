#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "types/p_dense_matrix/p_dense_matrix.h"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"

#include "types/utils/matrix_cast.h"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( trace_ambient, T, test_types)
{
    ambient::model >> dim(32,32), dim(32,32), dim(32,32);

    pMatrix pA(T::valuex,T::valuex);
    sMatrix sA(T::valuex,T::valuex);

    pA.set_init(ambient::random_i<typename T::dbl>);
    sA = maquis::traits::matrix_cast<sMatrix>(pA);

    typename T::dbl sa = trace(sA);
    typename T::dbl pa = trace(pA);

    maquis::cout << "Trace of sA " << sa << "; trace of pA " << pa << std::endl;
    Boost_check_close_adapter(sa,pa);
}


