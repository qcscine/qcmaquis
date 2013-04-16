#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( EXP, T, test_types)
{
    typedef alps::numeric::diagonal_matrix<typename T::value_type> sDiagMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<typename T::value_type> > pDiagMatrix;

    pDiagMatrix pA(T::valuex);
    sDiagMatrix sA((std::size_t)T::valuex);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sDiagMatrix>(pA);

    sA = exp(sA);
    ambient::numeric::exp_inplace(pA);

    BOOST_CHECK(pA==sA);
}

template<typename T>
struct random_helper{
    static void randomize(T& a){
        a = drand48();
    }

    static void randomize(std::complex<T>& a){
        a.real(drand48());
        a.imag(drand48());
    }
};


BOOST_AUTO_TEST_CASE_TEMPLATE( EXP_SCAL, T, test_types)
{
    typedef alps::numeric::diagonal_matrix<typename T::value_type> sDiagMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<typename T::value_type> > pDiagMatrix;

    typename T::value_type a; 
    random_helper<typename T::value_type>::randomize(a);

    pDiagMatrix pA(T::valuex);
    sDiagMatrix sA((std::size_t)T::valuex);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sDiagMatrix>(pA);

    sA = exp(sA*a);
    pA = expi(pA,a);

    BOOST_CHECK(pA==sA);
}
