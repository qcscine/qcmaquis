#include "utils/testing.hpp"

TEST_CASE( "Conjugate (inplace) is computed", "[conj_inplace]" )
{
    matrix<std::complex<double> >  A (TEST_M, TEST_N);
    matrix_<std::complex<double> > A_(TEST_M, TEST_N);

    generate(A);
    A_ = cast<matrix_<std::complex<double> > >(A);

    conj_inplace(A_);
    conj_inplace(A);

    REQUIRE((A == A_));
}

TEST_CASE( "Conjugate is computed", "[conj]" )
{
    matrix<std::complex<double> >  A (TEST_M, TEST_N);
    matrix_<std::complex<double> > A_(TEST_M, TEST_N);
    matrix<std::complex<double> >  B (TEST_M, TEST_N);
    matrix_<std::complex<double> > B_(TEST_M, TEST_N);

    generate(A);
    A_ = cast<matrix_<std::complex<double> > >(A);

    B_ = conj(A_);
    B  = conj(A);

    REQUIRE((A == A_));
    REQUIRE((B == B_));
}
