#include "utils/testing.hpp"

TEST_CASE( "Matrix inverse is computed", "[inverse]" )
{
    matrix<double> A(TEST_M, TEST_M);
    matrix<double> C(TEST_M, TEST_M);
    matrix<double> I(TEST_M, TEST_M);
    matrix<double> Id = matrix<double>::identity_matrix(TEST_M);
    
    generate(A);
    I = inverse(A);
    I *= A;
    
    REQUIRE((I == Id));
}
