#include "utils/testing.hpp"

template<class T>
struct test {
    static std::pair<matrix<T>, matrix_<T> > exponent(T alfa){
        matrix<T>  A (TEST_M, TEST_M), B (TEST_M, TEST_M);
        matrix_<T> A_(TEST_M, TEST_M), B_(TEST_M, TEST_M);
       
        generate_hermitian(A);
        A_ = cast<matrix_<T> >(A);
        
        B_ = exp_hermitian(A_, alfa);
        B  = exp_hermitian(A , alfa);

        return std::make_pair(B, B_);
    }
};

TEST_CASE( "Exponent (hermitian) is computed", "[exp_hermitian]" )
{
    double alfa;
    ambient::numeric::kernels::detail::randomize(alfa);
    const auto& pair = test<double>::exponent(alfa);
    REQUIRE((pair.first == pair.second));
}

TEST_CASE( "Exponent (hermitian, complex) is computed", "[exp_hermitian_complex]" )
{
    std::complex<double> alfa;
    ambient::numeric::kernels::detail::randomize(alfa);
    const auto& pair = test<std::complex<double> >::exponent(alfa);
    REQUIRE((pair.first == pair.second));
}

TEST_CASE( "Exponent is computed", "[exp]" )
{
    std::complex<double> a; 
    ambient::numeric::kernels::detail::randomize(a);
    
    matrix<std::complex<double> >  A (TEST_M, TEST_M), B (TEST_M, TEST_M);
    matrix_<std::complex<double> > A_(TEST_M, TEST_M), B_(TEST_M, TEST_M);
    
    generate(A);
    A_ = cast<matrix_<std::complex<double> > >(A);
    
    B_ = exp(A_, a);
    B  = exp(A , a);
    
    REQUIRE((B == B_));
}
