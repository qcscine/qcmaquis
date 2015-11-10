#include "utils/testing.hpp"

TEST_CASE( "Square root of diagonal matrix is computed", "[sqrt]" )
{
    diagonal<double> A(TEST_M, TEST_M);
    diagonal<double> S(TEST_M, TEST_M);
    diagonal<double> C(TEST_M, TEST_M);

    generate(A);
    S = sqrt(A);
    gemm(S, S, C);

    REQUIRE((A == C));
}
