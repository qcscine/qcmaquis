#include "utils/testing.hpp"

TEST_CASE( "Identity matrix is constructed", "[identity]" )
{
    matrix<double> A = matrix<double>::identity_matrix(TEST_M);
    double x = trace(A);
    double z = A(0,1);

    REQUIRE(x == TEST_M);
    REQUIRE(z == double());
}
