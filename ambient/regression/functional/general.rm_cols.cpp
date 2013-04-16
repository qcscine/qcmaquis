#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_FIRST_COL, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

//    sA.remove_cols(T::null,1);
//    pA.remove_cols(T::null,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_LAST_COL, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

//    sA.remove_cols(T::valuey-1,1);
//    pA.remove_cols(T::valuey-1,1);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_ONE_COL, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    int col =  Rd.IntRd()%(sA.num_cols());  
    
//    sA.remove_cols(col,1);
//    pA.remove_cols(col,1);

    BOOST_CHECK(pA==sA);
}
/*
BOOST_AUTO_TEST_CASE_TEMPLATE( REMOVE_SEVERAL_COLS, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);

    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
    int col =  Rd.IntRd()%(sA.num_cols()); 
    int numcols = T::valuey - col -1;   

    sA.remove_cols(col,numcols);
    pA.remove_cols(col,numcols);

    BOOST_CHECK(pA==sA);
}
*/
