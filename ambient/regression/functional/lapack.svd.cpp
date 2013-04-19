#include "params.hpp"

BOOST_AUTO_TEST_CASE_TEMPLATE( SVD_COMPARISON, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;
    typedef alps::numeric::diagonal_matrix<double> sDiagMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<double> > pDiagMatrix;

    pMatrix pA(T::valuex,T::valuey);
    pMatrix pU(T::valuex,T::valuey);
    pMatrix pV(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sU(T::valuex,T::valuey);
    sMatrix sV(T::valuex,T::valuey);

    pDiagMatrix pS(T::valuex);
    sDiagMatrix sS((std::size_t)T::valuex); 

    generate(pA);
    sA = cast<sMatrix>(pA);

    svd(pA,pU,pV,pS);
    svd(sA,sU,sV,sS);

    printf("Done %d x %d\n", (int)T::valuex, (int)T::valuey);
    BOOST_CHECK(sS == pS);
}
