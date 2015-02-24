#include "utils/testing.hpp"

TEST_CASE( "Conjugate is computed", "[conj]" )
{
    matrix<std::complex<double> >  A (TEST_M, TEST_N);
    matrix<std::complex<double> >  A_(TEST_M, TEST_N);

    generate(A);
    A_ = conj(A);
    conj_inplace(A_);

    REQUIRE((A == A_));
}
