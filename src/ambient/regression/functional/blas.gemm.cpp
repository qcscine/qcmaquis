#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( GEMM_NORMAL, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuex);
    pMatrix pB(T::valuex,T::valuex);
    pMatrix pC(T::valuex,T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sB(T::valuex,T::valuex);
    sMatrix sC(T::valuex,T::valuex);

    generate(pA);
    generate(pB);

    sA = matrix_cast<sMatrix>(pA);
    sB = matrix_cast<sMatrix>(pB);

    ambient::numeric::gemm(pA,pB,pC);
    ambient::sync();
    gemm(sA,sB,sC);
    BOOST_CHECK(pC==sC);
}

