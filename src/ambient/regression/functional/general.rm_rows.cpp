#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_FIRST_ROWS, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = matrix_cast<sMatrix>(pA);

//    sA.remove_rows(T::null,1);
//    pA.remove_rows(T::null,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_LAST_ROWS, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = matrix_cast<sMatrix>(pA);

//    sA.remove_rows(T::valuex-1,1);
//    pA.remove_rows(T::valuex-1,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_ROWS, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = matrix_cast<sMatrix>(pA);

    int row =  Rd.IntRd()%(sA.num_rows());  
//    sA.remove_rows(row,1);
//    pA.remove_rows(row,1);

    BOOST_CHECK(pA==sA);
}
/*
BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_SEVERAL_ROWS, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = matrix_cast<sMatrix>(pA);

    int row =  Rd.IntRd()%(T::valuex-1);  
    int numrows = (int)(T::valuex-1 - row)/2;   
    sA.remove_rows(row,numrows);
    pA.remove_rows(row,numrows);

    BOOST_CHECK(pA==sA);
}*/
