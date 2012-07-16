
#define BOOST_TEST_MODULE maquis::types::bench_qr
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/test/floating_point_comparison.hpp>

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

BOOST_AUTO_TEST_CASE_TEMPLATE(Bench_QRe, T, test_types)
{
    typedef maquis::types::dense_matrix<typename T::value_type> Matrix;

    Matrix A(T::valuex,T::valuey);
//    Matrix Q(T::valuex,T::valuey);
//    Matrix Qp(T::valuex,T::valuey);
    Matrix R(T::valuex,T::valuey);
//    Matrix RES(T::valuex,T::valuey);
    
    maquis::types::generate(A,drand48);
    Matrix D(A);
   
    std::string name("TimerQRe");
    name += "_" + boost::lexical_cast<std::string>(static_cast<int>(T::valuex)) + "x" + boost::lexical_cast<std::string>(static_cast<int>(T::valuex)); //need cast else crash compilation, ascii art
    int m(T::valuex);
    int info;
    int lwork(-1);
    double* work;
    double* work2;
    double wkopt;
    double tau[T::valuex];
    TimerOMP t(name);
    t.begin();
//    std::cout << A << std::endl;
    dgeqrf_(&m,&m, &A(0,0),&m,tau,&wkopt,&lwork,&info);
    work = new double[(int)wkopt];
    lwork = (int)wkopt;
    dgeqrf_(&m,&m, &A(0,0),&m,tau,work,&lwork,&info);

    for(int i=0; i < m ; ++i)
        for(int j=0; j <= i ; ++j)
            R(j,i) = A(j,i); 

    lwork = -1;
    dorgqr_(&m,&m,&m, &A(0,0),&m,tau,&wkopt,&lwork,&info);
    work2 = new double[(int)wkopt];
    lwork = (int)wkopt;
    dorgqr_(&m,&m,&m, &A(0,0),&m,tau,work2,&lwork,&info);
/*
    Qp = maquis::types::transpose(A);
    RES = A * Qp;
    std::cout << RES << std::endl;

    Q = A;
    Q =Q*R;
*/
    t.end();
/*
    for(int i=0; i < m ; ++i)
        for(int j=0; j < m ; ++j)
    BOOST_CHECK_CLOSE( Q(i,j), D(i,j), 0.0001 );
*/
    delete[] work;
    delete[] work2;

}
