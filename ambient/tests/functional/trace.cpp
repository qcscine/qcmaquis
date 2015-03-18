#include "utils/testing.hpp"

TEST_CASE( "Matrix trace is calculated", "[trace]" )
{
    matrix<double>  A(TEST_M,TEST_M);
    generate(A);
    double t = trace(A);

    double t_ = 0;
    for(int i = 0; i < TEST_M; i++) t_ += A(i,i);

    REQUIRE_CLOSE(t, t_);
}
