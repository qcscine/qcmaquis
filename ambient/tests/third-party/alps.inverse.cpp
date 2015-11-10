#include "utils/testing.hpp"

TEST_CASE( "Matrix inverse is computed", "[inverse]" )
{
    matrix<double>  A(TEST_M, TEST_M);
    matrix<double>  C(TEST_M, TEST_M);
    matrix<double>  I(TEST_M, TEST_M);
    matrix_<double> A_(TEST_M, TEST_M);
    matrix_<double> C_(TEST_M, TEST_M);
    matrix_<double> I_(TEST_M, TEST_M);
    
    generate(A);
    A_ = cast<matrix_<double> >(A);

    I  = inverse(A);
    I_ = inverse(A_);
    
    REQUIRE((I == I_));
}
