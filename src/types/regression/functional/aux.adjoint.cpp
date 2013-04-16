#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE(ADJOINT_INPLACE, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex, T::valuey);
    sMatrix sA(T::valuex, T::valuey);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    adjoint_inplace(sA);
    adjoint_inplace(pA);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ADJOINT, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex, T::valuey);
    pMatrix pB(T::valuex, T::valuey);
    sMatrix sA(T::valuex, T::valuey);
    sMatrix sB(T::valuex, T::valuey);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    sB = adjoint(sA);
    pB = adjoint(pA);

    BOOST_CHECK(pA==sA);
    BOOST_CHECK(pB==sB);
}
