#include "params.hpp"
#include "ambient/numeric/traits.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( GEMM_DIAGONAL, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;
    typedef alps::numeric::diagonal_matrix<typename T::value_type> sDiagMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<typename T::value_type> > pDiagMatrix;

    pMatrix pA(T::valuex,T::valuex);
    pMatrix pB(T::valuex,T::valuex);
    pDiagMatrix pC(T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sB(T::valuex,T::valuex);
    sDiagMatrix sC((std::size_t)T::valuex);

    generate(pB);
    generate(pC);

    sB = cast<sMatrix>(pB);
    sC = cast<sDiagMatrix>(pC);

    ambient::numeric::gemm(pC,pB,pA);
    gemm(sC,sB,sA);

    BOOST_CHECK(pA==sA);

    ambient::numeric::gemm(pB,pC,pA);
    gemm(sB,sC,sA);
    BOOST_CHECK(pA==sA);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( GEMM_DIAGONAL_COMPLEX_DOUBLE, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type > > pMatrix;

    typedef alps::numeric::diagonal_matrix< typename ambient::numeric::traits::real_type<T>::type  > sDiagMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<typename ambient::numeric::traits::real_type<T>::type > > pDiagMatrix;

    pMatrix pA(T::valuex,T::valuex);
    pMatrix pB(T::valuex,T::valuex);
    pDiagMatrix pC(T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sB(T::valuex,T::valuex);
    sDiagMatrix sC((std::size_t)T::valuex);

    generate(pB);
    generate(pC);

    sB = cast<sMatrix>(pB);
    sC = cast<sDiagMatrix>(pC);

    ambient::numeric::gemm(pC,pB,pA);
    gemm(sC,sB,sA);

    BOOST_CHECK(pA==sA);

    ambient::numeric::gemm(pB,pC,pA);
    gemm(sB,sC,sA);
    BOOST_CHECK(pA==sA);
}
