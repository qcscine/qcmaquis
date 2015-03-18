#include "utils/testing.hpp"
#include "ambient/container/numeric/traits.hpp"

TEST_CASE( "Matrix multiplication (with diagonal rhs) is computed", "[gemm_diagonal]" )
{
    matrix<double> A(TEST_M,TEST_M);
    matrix<double> B(TEST_M,TEST_M);
    diagonal<double> C(TEST_M,TEST_M);

    matrix_<double> A_(TEST_M,TEST_M);
    matrix_<double> B_(TEST_M,TEST_M);
    diagonal_<double> C_((std::size_t)TEST_M);

    generate(B);
    generate(C);

    B_ = cast<matrix_<double> >(B);
    C_ = cast<diagonal_<double> >(C);

    gemm(C,  B,  A);
    gemm(C_, B_, A_);

    REQUIRE((A == A_));

    gemm(B,  C,  A);
    gemm(B_, C_, A_);

    REQUIRE((A == A_));
}

TEST_CASE( "Matrix (complex) multiplication (with diagonal rhs) is computed", "[gemm_diagonal_complex]" )
{
    matrix<std::complex<double> > A(TEST_M,TEST_M);
    matrix<std::complex<double> > B(TEST_M,TEST_M);
    diagonal<double> C(TEST_M,TEST_M);

    matrix_<std::complex<double> > A_(TEST_M,TEST_M);
    matrix_<std::complex<double> > B_(TEST_M,TEST_M);
    diagonal_<double> C_((std::size_t)TEST_M);

    generate(B);
    generate(C);

    B_ = cast<matrix_<std::complex<double> > >(B);
    C_ = cast<diagonal_<double> >(C);

    gemm(C,  B,  A);
    gemm(C_, B_, A_);

    REQUIRE((A == A_));

    gemm(B,  C,  A);
    gemm(B_, C_, A_);

    REQUIRE((A == A_));
}
