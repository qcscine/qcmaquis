#include "utils/testing.hpp"

TEST_CASE( "Casting tiles to alps::matrix", "[cast_to_alps]" )
{
    matrix<double>  A (TEST_M,TEST_N);
    matrix_<double> A_(TEST_M,TEST_N);

    generate(A);
    A_ = cast<matrix_<double> >(A);

    REQUIRE((A == A_));
}

TEST_CASE( "Casting alps::matrix to tiles", "[cast_from_alps]" )
{
    matrix<double>  A (TEST_M,TEST_N);
    matrix_<double> A_(TEST_M,TEST_N);

    generate(A_, ambient::utils::Rd);
    A = cast<matrix<double> >(A_);

    REQUIRE((A_ == A));
}

TEST_CASE( "Casting tiles to alps::diagonal", "[cast_to_alps_diag]" )
{
    diagonal<double>  A (TEST_M,TEST_M);
    diagonal_<double> A_((std::size_t)TEST_M);

    generate(A);
    A_ = cast<diagonal_<double> >(A);

    REQUIRE((A == A_));
}

TEST_CASE( "Casting alps::diagonal to tiles", "[cast_from_alps_diag]" )
{
    diagonal<double>  H (TEST_M,TEST_M);
    diagonal<double>  A (TEST_M,TEST_M);
    diagonal_<double> A_((std::size_t)TEST_M);
   
    generate(H);
    A_ = cast<diagonal_<double> >(H);
    A  = cast<diagonal<double> >(A_);

    REQUIRE((A_ == A));
}
