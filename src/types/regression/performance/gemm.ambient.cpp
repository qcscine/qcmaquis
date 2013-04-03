#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "ambient/utils/timings.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( test, T, test_types){

    ambient::timer gtime("total"); gtime.begin();

    typedef ambient::dim2 dim;
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    size_t x = get_input_x<T>();
    size_t y = get_input_y<T>();
    size_t nthreads = get_input_threads<T>();

    pMatrix pA(x, y);
    pMatrix pB(x, y);
    pMatrix pC(x, y);
    pMatrix pC_orig(x, y);

    sMatrix sA(x, y);
    sMatrix sB(x, y);

    generate(pA);
    generate(pB);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    sB = maquis::bindings::matrix_cast<sMatrix>(pB);
    ambient::sync();

    printf("ambient::gemm strassen...\n");
    ambient::numeric::gemm_strassen(AMBIENT_MOVE(pA), AMBIENT_MOVE(pB), AMBIENT_MOVE(pC)); 
    ambient::timer time("ambient::gemm_strassen"); time.begin();
    ambient::sync();
    time.end();

    printf("ambient::gemm...\n");
    ambient::numeric::gemm(pA, pB, pC_orig); 
    ambient::timer time_orig("ambient::gemm"); time_orig.begin();
    ambient::sync();
    time_orig.end();

    report(time_orig, GFlopsGemm, x, y, nthreads);
    report(time, GFlopsGemm, x, y, nthreads);

    gtime.end();
    std::cout << "Global time: " << gtime.get_time() << "\n";
    BOOST_CHECK(pC==pC_orig);
}

