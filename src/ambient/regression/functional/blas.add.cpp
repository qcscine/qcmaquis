#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( ADD, T, test_types)
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

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    sB = maquis::bindings::matrix_cast<sMatrix>(pB);
    sC = maquis::bindings::matrix_cast<sMatrix>(pC);

    sC = sA + sB;
    pC = pA + pB; 

    BOOST_CHECK(pC == sC);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( ADD_ASSIGN, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    pMatrix pB(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sB(T::valuex,T::valuey);

    generate(pA);
    generate(pB);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    sB = maquis::bindings::matrix_cast<sMatrix>(pB);
 
    sA += sB; 
    pA += pB; 

    BOOST_CHECK(pA == sA);
}

