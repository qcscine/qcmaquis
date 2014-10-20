#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_FIRST_COL, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = cast<sMatrix>(pA);

    remove_cols(sA,0,1);
    remove_cols(pA,0,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_LAST_COL, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = cast<sMatrix>(pA);

    remove_cols(sA,T::valuey-1,1);
    remove_cols(pA,T::valuey-1,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_ONE_COL, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = cast<sMatrix>(pA);
    int col =  Rd.IntRd()%(sA.num_cols());  
    
    remove_cols(sA,col,1);
    remove_cols(pA,col,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_SEVERAL_COLS, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = cast<sMatrix>(pA);
    int col =  Rd.IntRd()%(sA.num_cols()); 
    int numcols = T::valuey - col -1;   

    remove_cols(sA,col,numcols);
    remove_cols(pA,col,numcols);

    BOOST_CHECK(pA==sA);
}
