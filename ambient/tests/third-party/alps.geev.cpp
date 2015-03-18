#include "utils/testing.hpp"

#include <boost/numeric/bindings/lapack/driver/geev.hpp>

TEST_CASE( "Matrix eigenvalues and left/right eigenvectors are computed", "[geev]" )
{
    matrix<std::complex<double> > A(TEST_M, TEST_M);
    matrix<std::complex<double> > L(TEST_M, TEST_M);
    matrix<std::complex<double> > R(TEST_M, TEST_M);

    matrix_<std::complex<double> > A_(TEST_M, TEST_M);
    matrix_<std::complex<double> > L_(TEST_M, TEST_M);
    matrix_<std::complex<double> > R_(TEST_M, TEST_M);

    diagonal<std::complex<double> > S(TEST_M, TEST_M); 
    typename alps::numeric::associated_vector< matrix_< std::complex<double> > >::type Sv_(num_rows(A_));

    generate(A);
    A_ = cast<matrix_<std::complex<double> > >(A);

    geev(A, L, R, S);

    int info = boost::numeric::bindings::lapack::geev('N', 'V', A_, Sv_, L_, R_);
    REQUIRE(info == 0);

    typename alps::numeric::associated_diagonal_matrix< matrix_< std::complex<double> > >::type S_(Sv_);

    REQUIRE((S == S_));
    REQUIRE((R == R_));
}
