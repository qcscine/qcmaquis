#include "utils/testing.hpp"

TEST_CASE( "Matrix transpose (inplace) is computed", "[transpose_inplace]" )
{
    matrix<double>  A (TEST_M,TEST_N);
    matrix_<double> A_(TEST_M,TEST_N);

    generate(A);
    A_ = cast<matrix_<double> >(A);

    transpose_inplace(A); 
    transpose_inplace(A_); 

    REQUIRE((A == A_));
}

