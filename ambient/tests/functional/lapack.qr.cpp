#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( QR_COMPARISON, T, test_types)
{
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix A(T::valuex,T::valuex);
    pMatrix Q(T::valuex,T::valuex);
    pMatrix R(T::valuex,T::valuex);
    pMatrix C(T::valuex,T::valuex);

    generate(A);
    qr(A,Q,R);
    gemm(Q,R,C);

    BOOST_CHECK(C == A);
}
