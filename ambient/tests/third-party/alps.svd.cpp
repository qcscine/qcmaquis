#include "utils/testing.hpp"

TEST_CASE( "Matrix singular value decomposition is performed", "[svd]" )
{
    matrix<double> A(TEST_M,TEST_N);
    matrix<double> U(TEST_M,TEST_N);
    matrix<double> V(TEST_M,TEST_N);

    matrix_<double> A_(TEST_M,TEST_N);
    matrix_<double> U_(TEST_M,TEST_N);
    matrix_<double> V_(TEST_M,TEST_N);

    diagonal<double>  S(TEST_M,TEST_M);
    diagonal_<double> S_((size_t)TEST_M); 

    generate(A);
    A_ = cast<matrix_<double> >(A);

    svd(A,  U,  V,  S);
    svd(A_, U_, V_, S_);

    REQUIRE((S_ == S));
}
