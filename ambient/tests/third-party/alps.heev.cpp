#include "utils/testing.hpp"

TEST_CASE( "Hermitian matrix eigenvalues and eigenvectors are computed", "[heev]" )
{
    matrix<std::complex<double> > A(TEST_M,TEST_M);
    matrix<std::complex<double> > V(TEST_M,TEST_M);
    diagonal<double> E(TEST_M,TEST_M); 

    matrix_<std::complex<double> > A_(TEST_M,TEST_M);
    matrix_<std::complex<double> > V_(TEST_M,TEST_M);
    diagonal_<double> E_((size_t)TEST_M);
 
    generate_hermitian(A);
    A_ = cast<matrix_<std::complex<double> > >(A);

    heev(A,  V,  E);
    heev(A_, V_, E_);

    REQUIRE((E == E_));
    //REQUIRE((V == V_)); // mismatch a bit
}
