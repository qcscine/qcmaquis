#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( HEEV_COMPARISON_VALUE, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<double> > pDiagMatrix;
    typedef alps::numeric::diagonal_matrix<double> sDiagMatrix;

    pMatrix pA(T::valuex,T::valuex);
    pMatrix pV(T::valuex,T::valuex);

    sMatrix sA(T::valuex,T::valuex);
    sMatrix sV(T::valuex,T::valuex);

    pDiagMatrix pE(T::valuex); 
    sDiagMatrix sE((std::size_t)T::valuex);
 
    generate(pA);
    sA = cast<sMatrix>(pA);

    heev(pA,pV,pE);
    heev(sA,sV,sE);

   // BOOST_CHECK(pE == sE);
    BOOST_CHECK(pV == sV);
}
