#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( CONJ_INPLACE, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex, T::valuey);
    sMatrix sA(T::valuex, T::valuey);

    if(T::valuex == T::valuey)
        generate_hermitian(pA);
    else
        generate(pA);

    sA = cast<sMatrix>(pA);

    bool bsA = is_hermitian(sA);
    bool bpA = is_hermitian(pA);

    BOOST_CHECK(bpA==bsA);
}

