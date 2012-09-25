#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <omp.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "ambient/numeric/matrix.hpp"
#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( test, T, test_types){
    typedef ambient::dim2 dim;
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    size_t x = get_input_x<T>();
    size_t y = get_input_y<T>();
    size_t nthreads = get_input_threads<T>();

    omp_set_num_threads(nthreads);
    pMatrix pA(x, y);
    pMatrix pL(x, y);
    pMatrix pQ(x, y);

    sMatrix sA(x, y);
    sMatrix sL(x, y);
    sMatrix sQ(x, y);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    ambient::sync();
    __a_timer time("ambient");
    time.begin();
    lq(sA, sL, sQ); 
    time.end();

    report(time, GFlopsGemm, x, y, nthreads);
}
