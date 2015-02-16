#include "utils/testing.hpp"

TEST_CASE( "Matrix multiplication is computed", "[gemm_normal]" )
{
    matrix<double> A(TEST_M,TEST_M);
    matrix<double> B(TEST_M,TEST_M);
    matrix<double> C(TEST_M,TEST_M);

    matrix_<double> A_(TEST_M,TEST_M);
    matrix_<double> B_(TEST_M,TEST_M);
    matrix_<double> C_(TEST_M,TEST_M);

    generate(A);
    generate(B);

    A_ = cast<matrix_<double> >(A);
    B_ = cast<matrix_<double> >(B);

    gemm(A,  B,  C);
    gemm(A_, B_, C_);

    REQUIRE((C == C_));
}

