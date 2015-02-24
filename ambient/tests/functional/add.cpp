#include "utils/testing.hpp"

TEST_CASE( "Matrix sum is computed", "[add]" )
{
    matrix<double> A(TEST_M,TEST_N);
    matrix<double> B(TEST_M,TEST_N);
    matrix<double> C(TEST_M,TEST_N);

    generate(A);
    generate(B);

    C = A + B;
    C -= B;

    REQUIRE((C == A));
}
