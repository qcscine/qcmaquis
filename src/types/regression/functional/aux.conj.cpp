
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/random.hpp>

#include "alps/numeric/matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( CONJ_INPLACE, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex, T::valuey);
    sMatrix sA(T::valuex, T::valuey);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    conj_inplace(sA);
    conj_inplace(pA);

    BOOST_CHECK(pA==sA);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CONJ, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex, T::valuey);
    pMatrix pB(T::valuex, T::valuey);
    sMatrix sA(T::valuex, T::valuey);
    sMatrix sB(T::valuex, T::valuey);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    sB = conj(sA);
    pB = conj(pA);

    BOOST_CHECK(pA==sA);
    BOOST_CHECK(pB==sB);
}
