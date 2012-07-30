#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( remove_first_col, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    pA.fill_random();

    sA = maquis::bindings::matrix_cast<sMatrix>(pA); // playout is inside the cast

    sA.remove_cols(T::null,1);
    pA.remove_cols(T::null,1);

    BOOST_CHECK(pA==sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

BOOST_AUTO_TEST_CASE_TEMPLATE( remove_last_col, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    pA.fill_random();

    sA = maquis::bindings::matrix_cast<sMatrix>(pA); // playout is inside the cast

    sA.remove_cols(T::valuey-1,1);
    pA.remove_cols(T::valuey-1,1);

    BOOST_CHECK(pA==sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}

BOOST_AUTO_TEST_CASE_TEMPLATE( remove_one_col, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    pA.fill_random();

    sA = maquis::bindings::matrix_cast<sMatrix>(pA); // playout is inside the cast
    int col =  Rd.IntRd()%(sA.num_cols());  
    
    sA.remove_cols(col,1);
    pA.remove_cols(col,1);

    BOOST_CHECK(pA==sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}
/*
BOOST_AUTO_TEST_CASE_TEMPLATE( remove_several_cols, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    pA.fill_random();

    sA = maquis::bindings::matrix_cast<sMatrix>(pA); // playout is inside the cast
    int col =  Rd.IntRd()%(sA.num_cols()); 
    int numcols = T::valuey - col -1;   

    sA.remove_cols(col,numcols);
    pA.remove_cols(col,numcols);

    BOOST_CHECK(pA==sA); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}
*/
