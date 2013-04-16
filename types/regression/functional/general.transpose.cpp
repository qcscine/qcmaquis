#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( TRANSPOSE, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    transpose_inplace(pA); 
    transpose_inplace(sA); 

    BOOST_CHECK(pA==sA);
}

