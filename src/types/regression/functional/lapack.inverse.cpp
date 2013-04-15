
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

#include <boost/numeric/bindings/lapack/driver/geev.hpp>

// this function call geqrf and geqri so lapack test
BOOST_AUTO_TEST_CASE_TEMPLATE( INVERSE_COMPARISON, T, test_types)
{
   // test_complex_only<typename T::value_type>::test(T::valuex);
        typedef alps::numeric::matrix<typename T::value_type> sMatrix;
        typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;
       
        int size = T::valuex;
 
        pMatrix  pA(size, size);
        pMatrix  pC(size, size);
        pMatrix  pAinv(size, size);
        sMatrix  sA(size, size);
        sMatrix  sC(size, size);
        sMatrix  sAinv(size, size);
        
        generate(pA);
        sA = maquis::bindings::matrix_cast<sMatrix>(pA);
        pAinv = inverse(pA);
        sAinv = inverse(sA);
        
        BOOST_CHECK(pAinv == sAinv);
}
