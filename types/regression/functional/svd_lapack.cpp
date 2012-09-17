
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "ambient/numeric/matrix.hpp"
#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( addition, T, test_types)
{
    pMatrix pA(T::valuex,T::valuey);
    pMatrix pU(T::valuex,T::valuey);
    pMatrix pV(T::valuex,T::valuey);

    sMatrix sA(T::valuex,T::valuey);
    sMatrix sU(T::valuex,T::valuey);
    sMatrix sV(T::valuex,T::valuey);

    typename alps::numeric::associated_real_diagonal_matrix<pMatrix>::type pS;
    typename alps::numeric::associated_real_diagonal_matrix<sMatrix>::type sS; 

    fill_random(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);
 
    svd(pA,pU,pV,pS);
    svd(sA,sU,sV,sS);
  
    BOOST_CHECK(sS == pS); // BOOST_CHECK_EQUAL necessitates == inside the class, here == is a free function 
}
