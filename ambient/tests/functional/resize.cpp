#include "utils/testing.hpp"

TEST_CASE( "Matrix is constructed and resized", "[resize]" )
{
    matrix<double> A = matrix<double>::identity_matrix(TEST_M);
    matrix<double> H = matrix<double>::identity_matrix(TEST_M/2);

    A.resize(TEST_M,TEST_M/2);

    REQUIRE((A.num_rows() == TEST_M));
    REQUIRE((A.num_cols() == TEST_M/2));

    A.resize(TEST_M/2,TEST_M/2);
    REQUIRE((A == H));
}
