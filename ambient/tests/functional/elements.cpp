#include "utils/testing.hpp"

TEST_CASE( "Element access is performed", "[element_access]" )
{
    size_t accessx = TEST_M-1;
    size_t accessy = TEST_N-1;

    matrix<double> A(TEST_M,TEST_N);
    generate(A);
    A (accessx,accessy) = 3;

    REQUIRE((A(accessx,accessy) == 3));
}
