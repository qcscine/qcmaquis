#include "utils/testing.hpp"

TEST_CASE( "Matrix adjoint is computed (inplace)", "[adjoint_inplace]" )
{
    matrix<std::complex<double> >  A (TEST_M, TEST_N);
    matrix_<std::complex<double> > A_(TEST_M, TEST_N);

    generate(A);
    A_ = cast<matrix_<std::complex<double> > >(A);

    adjoint_inplace(A_);
    adjoint_inplace(A);

    REQUIRE((A == A_));
}

TEST_CASE( "Matrix adjoint is computed", "[adjoint]" )
{
    matrix<std::complex<double> >  A (TEST_M, TEST_N);
    matrix_<std::complex<double> > A_(TEST_M, TEST_N);
    matrix<std::complex<double> >  B (TEST_M, TEST_N);
    matrix_<std::complex<double> > B_(TEST_M, TEST_N);

    generate(A);
    A_ = cast<matrix_<std::complex<double> > >(A);

    B_ = adjoint(A_);
    B  = adjoint(A);

    REQUIRE((A==A_));
    REQUIRE((B==B_));
}
