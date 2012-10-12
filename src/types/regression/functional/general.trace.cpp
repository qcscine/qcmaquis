#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "ambient/numeric/matrix.hpp"
#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE( TRACE, T, test_types)
{
    pMatrix pA(T::valuex,T::valuex);
    sMatrix sA(T::valuex,T::valuex);

    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    typename T::dbl sa = trace(sA);
    typename T::dbl pa = trace(pA);

    maquis::cout << "Trace of sA " << sa << "; trace of pA " << pa << std::endl;
    Boost_check_close_adapter(sa,pa);
}


