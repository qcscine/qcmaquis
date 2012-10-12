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

BOOST_AUTO_TEST_CASE_TEMPLATE( BIDIAGONALIZATION, T, test_types)
{
    pMatrix A(T::valuex,T::valuey);
    pMatrix U(T::valuex,T::valuey);
    pMatrix V(T::valuex,T::valuey);
    pMatrix S(T::valuex,T::valuey);

    generate(A);
    ambient::numeric::gebrd(A,U,S,V);

    pMatrix G1(num_rows(U),num_cols(S));
    pMatrix G2(num_rows(U),num_cols(V));

    ambient::numeric::gemm(U, S, G1);
    ambient::numeric::gemm(G1, V, G2);

    std::cout << "done: " << num_rows(A) << " " << num_cols(A) << "\n";
    BOOST_CHECK(G2 == A);
}
