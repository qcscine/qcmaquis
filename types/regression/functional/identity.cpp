

#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/random.hpp>


#include "types/p_dense_matrix/p_dense_matrix.h"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"

#include "types/utils/bindings.hpp"

#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( identity, T, test_types)
{
    ambient::model >> dim(32,32), dim(32,32), dim(32,32);
    std::size_t accessx(T::valuex-1), accessy(T::valuey-1);
    typename T::dbl x;
    typename T::dbl y;
    // check if we are writing inside the matrix
    BOOST_STATIC_ASSERT(T::valuex-1>0);
    BOOST_STATIC_ASSERT(T::valuex-1>0);

    pMatrix pA(T::valuex,T::valuey);
    sMatrix sA(T::valuex,T::valuey);

    maquis::types::algorithms::generate(sA,Rd); // Rd is rand generator static variable inside utilities
    pA = maquis::traits::matrix_cast<pMatrix>(sA); // playout is inside the cast

    pA(accessx,accessy) = 3;
    ambient::playout();
    sA(accessx,accessy) = 3;

    x =  pA(accessx,accessy) ;
    y =  sA(accessx,accessy) ;

    Boost_check_close_adapter(x,y);
    BOOST_CHECK(pA==sA); // we also check everything to verify we do not corrupt the memory
}

BOOST_AUTO_TEST_CASE_TEMPLATE( read_access, T, test_types)
{
    ambient::model >> dim(32,32), dim(32,32), dim(32,32);
    typename T::dbl x,y;
    pMatrix pA = pMatrix::identity_matrix(T::valuex);
    sMatrix sA = sMatrix::identity_matrix(T::valuex);
   
    x = trace(pA);
    y = trace(sA);
   
    ambient::playout();
    Boost_check_close_adapter(x,y);
}

