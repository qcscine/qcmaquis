
#define BOOST_TEST_MODULE maquis::types::bench_qr
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

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


// Define a base random number generator and initialize it with a seed.
boost::random::mt19937 rng(3); 
// Define distribution U[0,1) [double values]
boost::random::uniform_real_distribution<> dist(0,1);
// Define a random variate generator using our base generator and distribution
boost::variate_generator<boost::random::mt19937&, boost::random::uniform_real_distribution<> > uniDblGen(rng, dist);


BOOST_AUTO_TEST_CASE_TEMPLATE(QR_test, T, test_types)
{
    typedef maquis::types::dense_matrix<typename T::value_type> Matrix;

    Matrix M(T::valuex,T::valuey);
    Matrix Q;
    Matrix R;
    Matrix A;



    maquis::types::generate(M,uniDblGen);
    std::cout << M << std::endl;
    std::string name("TimerQR");
    name += "_" + boost::lexical_cast<std::string>(static_cast<int>(T::valuex)) + "x" + boost::lexical_cast<std::string>(static_cast<int>(T::valuey)); //need cast else crash compilation, ascii art
    TimerOMP t(name);
    t.begin();
    maquis::types::qr(M,Q,R);
    t.end();

    std::cout << num_rows(M) << " " << num_cols(M) << std::endl;
    std::cout << num_rows(Q) << " " << num_cols(Q) << std::endl;
    std::cout << num_rows(R) << " " << num_cols(R) << std::endl;
 // first test we check Q
    Matrix D(Q*R);

    for(int i(0); i< num_rows(M); ++i) 
        for(int j(0); j< num_cols(M); ++j) 
           BOOST_CHECK_CLOSE(M(i,j),D(i,j),1e-6); 



}

BOOST_AUTO_TEST_CASE_TEMPLATE(Q_ID_test, T, test_types)
{
    typedef maquis::types::dense_matrix<typename T::value_type> Matrix;

    Matrix M(T::valuex,T::valuey);
    Matrix Q;
    Matrix R;
    Matrix A;

    maquis::types::generate(M,uniDblGen);
    std::string name("TimerQR");
    name += "_" + boost::lexical_cast<std::string>(static_cast<int>(T::valuex)) + "x" + boost::lexical_cast<std::string>(static_cast<int>(T::valuey)); //need cast else crash compilation, ascii art
    TimerOMP t(name);
    t.begin();
    maquis::types::qr(M,Q,R);
    t.end();

 // second test we check Q is an Id matrix, cautions we implemented the thin QR so onli Qt*Q is equal to one
    Matrix D(transpose(Q)*Q);

    for(int i(0); i< std::min(num_rows(D),num_cols(D)); ++i) 
        BOOST_CHECK_CLOSE(1.0,D(i,i),1e-6); 
}

