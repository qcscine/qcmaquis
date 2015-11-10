#include "utils/testing.hpp"

TEST_CASE( "Matrix is Hermitian", "[is_hermitian]" )
{
    matrix<double> A(TEST_M, TEST_M);
    matrix<std::complex<double> > Ac(TEST_M, TEST_M);

    generate_hermitian(A);
    generate_hermitian(Ac);

    REQUIRE((is_hermitian(A)  == true));
    REQUIRE((is_hermitian(Ac) == true));
}

