#include "utils/testing.hpp"

TEST_CASE( "Matrix transpose (inplace) is computed", "[transpose_inplace]" )
{
    matrix<double> A (TEST_M,TEST_N);
    matrix<double> A_(TEST_M,TEST_N);

    generate(A);
    A_ = transpose(A); 
    transpose_inplace(A_);

    REQUIRE((A == A_));
}
