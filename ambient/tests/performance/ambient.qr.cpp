#include "utils/testing.hpp"

TEST_CASE( "Matrix QR factorization performance measured", "[ambient::qr]" )
{
    measurement params;
    size_t x = params.num_cols();
    size_t y = params.num_rows();

    matrix<double> A(x, y);
    matrix<double> Q(x, y);
    matrix<double> R(x, y);

    generate(A);
    ambient::sync();

    qr(A,  Q,  R); 

    measurement::timer time("qr"); time.begin();
    ambient::sync();
    time.end();

    params.report(gflops::gemm, time.get_time());
}

