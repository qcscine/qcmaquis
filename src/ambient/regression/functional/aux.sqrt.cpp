#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( SQRT, T, test_types)
{
    typedef alps::numeric::diagonal_matrix<typename T::value_type> sDiagMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<typename T::value_type> > pDiagMatrix;

    pDiagMatrix pA(T::valuex);
    sDiagMatrix sA((std::size_t)T::valuex);

    generate(pA);
    sA = matrix_cast<sDiagMatrix>(pA);

    sA = sqrt(sA);
    ambient::numeric::sqrt_inplace(pA);

    BOOST_CHECK(pA==sA);
}
