#include "utils/testing.hpp"

TEST_CASE( "Matrix adjoint is computed", "[adjoint]" )
{
    matrix<std::complex<double> > A (TEST_M, TEST_N);
    matrix<std::complex<double> > A_(TEST_N, TEST_M);

    generate(A);
    A_ = adjoint(A);
    adjoint_inplace(A_);

    REQUIRE((A == A_));
}
