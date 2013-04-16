#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( CAST_P2S, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);
    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( CAST_S2P, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);
    generate(sA,Rd); // Rd is rand generator static variable inside utilities
    pA = maquis::bindings::matrix_cast<pMatrix>(sA);
    BOOST_CHECK(sA==pA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( CAST_C2D, T, test_types) //complex to double only for TE into exp function heev
{
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<double> > pDiagMatrixDouble;
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix< std::complex<double> > > pDiagMatrixComplex;

    pDiagMatrixDouble pAd(T::valuex,T::valuey);
    pDiagMatrixDouble pBd(T::valuex,T::valuey);
    pDiagMatrixComplex pAc(T::valuex,T::valuey);

    generate(pAd,Rd); // Rd is rand generator static variable inside utilities
    pAc = maquis::bindings::matrix_cast<pDiagMatrixComplex,pDiagMatrixDouble>(pAd);
    pBd = maquis::bindings::matrix_cast<pDiagMatrixDouble,pDiagMatrixComplex>(pAc);
    BOOST_CHECK(pAd==pBd);
}


/*
BOOST_AUTO_TEST_CASE_TEMPLATE( CAST_P2S_DIAG, T, test_types)
{
    pDiagMatrix pA(T::valuex);
    sDiagMatrix sA((std::size_t)T::valuex);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sDiagMatrix>(pA);
    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( CAST_S2P_DIAG, T, test_types)
{
    pDiagMatrix pA(T::valuex);
    sDiagMatrix sA((std::size_t)T::valuex);
   
    sA.generate(Rd); // Rd is rand generator static variable inside utilities
    pA = maquis::bindings::matrix_cast<pDiagMatrix>(sA);
    BOOST_CHECK(sA==pA);
}
*/
