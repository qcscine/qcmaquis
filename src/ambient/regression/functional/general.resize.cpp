#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( RESIZE, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    generate(pA);
    sA = matrix_cast<sMatrix>(pA);

    sA.resize(4,2);
    pA.resize(4,2);

    BOOST_CHECK(pA == sA);
}
