#include "utils/testing.hpp"

TEST_CASE( "Matrix LQ factorization performance measured", "[lapack::lq]" )
{
    measurement params;

    size_t x = params.num_cols();
    size_t y = params.num_rows();

    matrix<double>  A (x, y);
    matrix_<double> A_(x, y);
    matrix_<double> L_(x, y);
    matrix_<double> Q_(x, y);

    generate(A);
    A_ = cast<matrix_<double> >(A);
    ambient::sync();

    measurement::timer time("lq"); time.begin();
    lq(A_, L_, Q_);
    time.end();

    params.report(gflops::gemm, time.get_time());
}
