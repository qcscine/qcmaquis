#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( QR_COMPARISON, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    pMatrix pQ(T::valuex,T::valuey);
    pMatrix pR(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sQ(T::valuex,T::valuey);
    sMatrix sR(T::valuex,T::valuey);


    generate(pA);
    sA = cast<sMatrix>(pA);

    qr(pA,pQ,pR);
    qr(sA,sQ,sR);

    BOOST_CHECK(sQ == pQ);
    BOOST_CHECK(sR == pR);
}
