#include "utils/testing.hpp"

TEST_CASE( "Matrix first column is removed", "[remove_first_col]" )
{
    matrix<double> A = matrix<double>::identity_matrix(TEST_M);
    remove_cols(A, 0, 1);

    REQUIRE((A(1,0) == (double)1));
}

TEST_CASE( "Matrix last column is removed", "[remove_last_col]" )
{
    matrix<double> A(TEST_M,TEST_N);
    generate(A);

    remove_cols(A, TEST_N-1, 1);
    ambient::sync();

    REQUIRE(( A.num_cols() == TEST_N-1 ));
}

TEST_CASE( "Matrix column is removed", "[remove_one_col]" )
{
    matrix<double> A = matrix<double>::identity_matrix(TEST_M);
    int col =  ambient::utils::Rd.IntRd() % A.num_cols();
    remove_cols(A, col, 1);

    REQUIRE((A(col+1,col) == (double)1));
}

TEST_CASE( "Some matrix columns are removed", "[remove_several_cols]" )
{
    matrix<double> A = matrix<double>::identity_matrix(TEST_M);
    int col =  ambient::utils::Rd.IntRd() % (TEST_M-1);
    int numcols = TEST_M - col - 1;
    remove_cols(A, col, numcols);

    REQUIRE((A(col+numcols,col) == (double)1));
}
