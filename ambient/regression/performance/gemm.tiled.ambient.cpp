#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( test, T, test_types){
    typedef ambient::dim2 dim;
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles< ambient::numeric::matrix<typename T::value_type> > pMatrix;

    size_t x = get_input_x<T>();
    size_t y = get_input_y<T>();
    size_t nthreads = get_input_threads<T>();

    pMatrix pA(x, y);
    pMatrix pB(x, y);
    pMatrix pC(x, y);

    sMatrix sA(x, y);
    sMatrix sB(x, y);

    generate(pA);
    generate(pB);

    sA = matrix_cast<sMatrix>(pA);
    sB = matrix_cast<sMatrix>(pB);
    ambient::sync();

    ambient::numeric::gemm(pA, pB, pC); 

    ambient::timer time("ambient");
    time.begin();
    ambient::sync();
    time.end();

    report(time, GFlopsGemm, x, y, nthreads);
}

