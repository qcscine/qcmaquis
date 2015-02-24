#include "utils/testing.hpp"

TEST_CASE( "Matrix QR factorization performance measured", "[lapack::qr]" )
{
    measurement params;

    size_t x = params.num_cols();
    size_t y = params.num_rows();

    matrix<double>  A (x, y);
    matrix_<double> A_(x, y);
    matrix_<double> Q_(x, y);
    matrix_<double> R_(x, y);

    generate(A);
    A_ = cast<matrix_<double> >(A);
    ambient::sync();

    measurement::timer time("qr"); time.begin();
    qr(A_, Q_, R_); 
    time.end();

    params.report(gflops::gemm, time.get_time());
}
