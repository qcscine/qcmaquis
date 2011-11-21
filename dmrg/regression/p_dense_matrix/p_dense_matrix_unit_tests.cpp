#include "utils/zout.hpp"
#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "types/p_dense_matrix/concept/matrix_interface.hpp"
#include "types/p_dense_matrix/concept/resizable_matrix_interface.hpp"

#define BOOST_TEST_MODULE p_dense_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>

#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/syevd.hpp>

#define M_SIZE 128
#define N_SIZE 256
using namespace blas;

typedef boost::mpl::list<double> test_types;
typedef ambient::dim2 dim;

struct caveats {
    caveats(){ ambient::init();     }
   ~caveats(){ ambient::finalize(); }
};

BOOST_GLOBAL_FIXTURE( caveats );
/*
BOOST_AUTO_TEST_CASE_TEMPLATE( identity_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 

    p_dense_matrix<T> A = p_dense_matrix<T>::identity_matrix(2);
    __ambient_wo_begin__
    A(1,0) = 3;
    __ambient_wo_end__
    std::cout << A;
    p_dense_matrix<T> B(A);
    B.resize(10,10);
    __ambient_wo_begin__
    B(0,1) = 2;
    B(5,0) = 26;
    __ambient_wo_end__
    std::cout << B;
}
BOOST_AUTO_TEST_CASE_TEMPLATE( scalar_norm_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 

    size_t length = 4;
    std::vector< p_dense_matrix<double> > data_;
    std::vector< p_dense_matrix<double> > ret;

    data_.reserve(length); // avoid calling copy constructors on PUSH_BACK (init will reset)
    ret.reserve(length); // avoid calling copy constructors on PUSH_BACK (init will reset)
    for(size_t i = 0 ; i < length; i++){
        p_dense_matrix<double> A(6,6);
        data_.push_back(A);
        generate(data_[i], double());
    }
    for(size_t i = 0 ; i < length; i++){
        ret.push_back(p_dense_matrix<double>(1,1));
        ambient::push(ambient::scalar_norm_l, ambient::scalar_norm_c, data_[i], ret[i]);
        printf("PARTIAL NORM: %.6f\n", ret[i](0,0));
    }
    for(size_t i = 1 ; i < length; i++){
        ambient::push(ambient::atomic_add_l, ambient::atomic_add_c, ret[0], ret[i]);
    }
    printf("TOTAL NORM: %.6f\n", ret[0](0,0));

}
*/
BOOST_AUTO_TEST_CASE_TEMPLATE( bela_test, T, test_types ) 
{ 
   typedef p_dense_matrix<T, ambient::MANUAL> MatrixT;
   ambient::layout >> dim(1,1), dim(1,1), dim(10,1);

   MatrixT a(2,2);
   a(0,0) = 1;

   std::vector<MatrixT> foo(3, a);
   for (int i = 0; i < 3; ++i)
       std::cout << foo[i](0,0) << std::endl;

   foo[1](0,0) = 2;
   for (int i = 0; i < 3; ++i)
       std::cout << foo[i](0,0) << std::endl;

   foo.push_back(MatrixT(3,3));
   foo[3](0,0) = 3;
   foo[0](0,0) = -1;
   for (int i = 0; i < 4; ++i)
       std::cout << foo[i](0,0) << std::endl;

}


