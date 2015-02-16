#include "utils/testing.hpp"

TEST_CASE( "Matrix first rows are removed", "[remove_first_rows]" )
{
    matrix<double> A = matrix<double>::identity_matrix(TEST_M);
    remove_rows(A, 0, 1);

    REQUIRE((A(0,1) == (double)1));
}

TEST_CASE( "Matrix last row is removed", "[remove_last_rows]" )
{
    matrix<double>  A(TEST_M,TEST_N);
    generate(A);

    remove_rows(A, TEST_M-1, 1);
    ambient::sync();

    REQUIRE((A.num_rows() == TEST_M-1));
}

TEST_CASE( "Matrix row is removed", "[remove_row]" )
{
    matrix<double> A = matrix<double>::identity_matrix(TEST_M);
    int row = ambient::utils::Rd.IntRd() % A.num_rows();
    remove_rows(A, row, 1);

    REQUIRE((A(row,row+1) == (double)1));
}

TEST_CASE( "Several matrix rows are removed", "[remove_several_rows]" )
{
    matrix<double> A = matrix<double>::identity_matrix(TEST_M);
    int row = ambient::utils::Rd.IntRd() % (TEST_M-1);
    int numrows = (int)(TEST_M - 1 - row)/2;
    remove_rows(A, row, numrows);

    REQUIRE((A(row,row+numrows) == (double)1));
}
