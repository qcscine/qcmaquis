#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_FIRST_ROWS, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = cast<sMatrix>(pA);

    remove_rows(sA,0,1);
    remove_rows(pA,0,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_LAST_ROWS, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = cast<sMatrix>(pA);

    remove_rows(sA,T::valuex-1,1);
    remove_rows(pA,T::valuex-1,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_ROWS, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = cast<sMatrix>(pA);

    int row =  Rd.IntRd()%(sA.num_rows());  
    remove_rows(sA,row,1);
    remove_rows(pA,row,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_SEVERAL_ROWS, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = cast<sMatrix>(pA);

    int row =  Rd.IntRd()%(T::valuex-1);  
    int numrows = (int)(T::valuex-1 - row)/2;   
    remove_rows(sA,row,numrows);
    remove_rows(pA,row,numrows);

    BOOST_CHECK(pA==sA);
}
