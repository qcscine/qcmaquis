#include "utils/testing.hpp"

TEST_CASE( "Matrix LQ factorization performance measured", "[ambient::lq]" )
{
    measurement params;
    size_t x = params.num_cols();
    size_t y = params.num_rows();

    matrix<double> A(x, y);
    matrix<double> Q(x, y);
    matrix<double> L(x, y);

    generate(A);
    ambient::sync();

    lq(A, L, Q); 

    measurement::timer time("lq"); time.begin();
    ambient::sync();
    time.end();

    params.report(gflops::gemm, time.get_time());
}

