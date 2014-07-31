#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( SUB, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    pMatrix pB(T::valuex,T::valuey);
    pMatrix pC(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sB(T::valuex,T::valuey);
    sMatrix sC(T::valuex,T::valuey);

    generate(pA);
    generate(pB);

    sA = cast<sMatrix>(pA);
    sB = cast<sMatrix>(pB);
    sC = cast<sMatrix>(pC);
 
    sC = sA - sB; 
    pC = pA - pB; 

    BOOST_CHECK(pC == sC);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SUB_ASSIGN, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    pMatrix pB(T::valuex,T::valuey);
    pMatrix pC(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sB(T::valuex,T::valuey);
    sMatrix sC(T::valuex,T::valuey);

    generate(pA);
    generate(pB);

    sA = cast<sMatrix>(pA);
    sB = cast<sMatrix>(pB);
    sC = cast<sMatrix>(pC);
 
    sA -= sB; 
    pA -= pB; 

    BOOST_CHECK(pA == sA);
}

