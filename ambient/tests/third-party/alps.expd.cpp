#include "utils/testing.hpp"

template<class T>
struct test {
    static std::pair<diagonal<T>, diagonal_<T> > exponent(){
        diagonal<T>  A (TEST_M, TEST_M);
        diagonal_<T> A_((std::size_t)TEST_M);
        generate(A); A_ = cast<diagonal_<T> >(A);
        
        A_ = exp(A_);
        exp_inplace(A);
        return std::make_pair(A, A_);
    }

    static std::pair<diagonal<T>, diagonal_<T> > exponent(T alfa){
        diagonal<T>  A (TEST_M, TEST_M);
        diagonal_<T> A_((std::size_t)TEST_M);
        
        generate(A);
        A_ = cast<diagonal_<T> >(A);
        
        A_ = exp(A_*alfa);
        A  = expi(A,alfa);
        return std::make_pair(A, A_);
    }
};

TEST_CASE( "Exponential is computed (diagonal)", "[exp_diagonal]" )
{
    const auto& pair = test<double>::exponent();
    REQUIRE((pair.first == pair.second));

    const auto& pairc = test<std::complex<double> >::exponent();
    REQUIRE((pair.first == pair.second));
}

TEST_CASE( "Exponential (with factor) is computed (diagonal)", "[exp_diagonal_scal]" )
{
    std::complex<double> alfa;
    ambient::numeric::kernels::detail::randomize(alfa);

    const auto& pair = test<double>::exponent(alfa.real());
    REQUIRE((pair.first == pair.second));

    const auto& pairc = test<std::complex<double> >::exponent(alfa);
    REQUIRE((pair.first == pair.second));
}
