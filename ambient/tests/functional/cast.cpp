#include "utils/testing.hpp"

TEST_CASE( "Casting complex matrix to real matrix", "[cast_complex_double]" )
{
    diagonal<double> A(TEST_M,TEST_M);
    diagonal<double> B(TEST_M,TEST_M);
    diagonal<std::complex<double> > Ac(TEST_M,TEST_M);

    generate(A);
    Ac = cast<diagonal<std::complex<double> >, diagonal<double> >(A);
    B = cast<diagonal<double>, diagonal<std::complex<double> > >(Ac);
    REQUIRE((A == B));
}
