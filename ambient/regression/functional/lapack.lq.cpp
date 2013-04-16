#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( LQ_COMPARISON, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    pMatrix pL(T::valuex,T::valuey);
    pMatrix pQ(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sL(T::valuex,T::valuey);
    sMatrix sQ(T::valuex,T::valuey);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
 
    lq(pA,pL,pQ);
    lq(sA,sL,sQ);

    BOOST_CHECK(sL == pL);
    BOOST_CHECK(sQ == pQ);
}
