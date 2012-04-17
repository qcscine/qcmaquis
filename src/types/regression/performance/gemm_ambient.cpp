#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>
#include <cmath>

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/algorithms.hpp"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/utils/matrix_cast.h"

#include "utils/timings.h"
#include "utilities.h"


BOOST_AUTO_TEST_CASE_TEMPLATE( test, T, test_types){
    typedef ambient::dim2 dim;
    typedef maquis::types::dense_matrix<typename T::value_type> sMatrix;
    typedef maquis::types::diagonal_matrix<typename T::value_type> sDiagMatrix;
    typedef maquis::types::p_dense_matrix<typename T::value_type> pMatrix;
    typedef maquis::types::p_diagonal_matrix<typename T::value_type> pDiagMatrix;

    size_t x = get_input_x<T>();
    size_t y = get_input_y<T>();
    size_t nthreads = get_input_threads<T>();

    ambient::model >> dim(256,256), dim(256,256), dim(256,256);
    ambient::set_num_threads(nthreads);

    pMatrix pA(x, y);
    pMatrix pB(x, y);
    pMatrix pC(x, y);

    sMatrix sA(x, y);
    sMatrix sB(x, y);

    pA.set_init(ambient::random_i<typename T::value_type>);
    pB.set_init(ambient::random_i<typename T::value_type>);

    sA = maquis::traits::matrix_cast<sMatrix>(pA); // playout is inside the cast
    sB = maquis::traits::matrix_cast<sMatrix>(pB); // playout is inside the cast
    ambient::playout();

    maquis::types::algorithms::gemm(pA, pB, pC); 

    TimerPTH time("ambient");
    time.begin();
    ambient::playout();
    time.end();

    report(time, GFlopsGemm, x, y, nthreads);
}

