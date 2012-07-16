
#define BOOST_TEST_MODULE maquis::types::bench_svd
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/lexical_cast.hpp>

#include <string>

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/algorithms.hpp"

#include "types/dense_matrix/vector.hpp"
#include "utils/timings.h"
#include "utilities.h"

BOOST_AUTO_TEST_CASE_TEMPLATE(Bench_SVD, T, test_types)
{
    typedef maquis::types::dense_matrix<typename T::value_type> Matrix;

    Matrix A(T::valuex,T::valuey);
    Matrix U(T::valuex,T::valuey);
    Matrix V(T::valuex,T::valuey);

    typename maquis::types::associated_diagonal_matrix<Matrix>::type S; 
    
    maquis::types::generate(A,drand48);
    std::string name("TimerSVD");
    name += "_" + boost::lexical_cast<std::string>(static_cast<int>(T::valuex)) + "x" + boost::lexical_cast<std::string>(static_cast<int>(T::valuex)); //need cast else crash compilation, ascii art
    TimerOMP t(name);
    t.begin();
    maquis::types::svd(A,U,V,S);
    t.end();
}
