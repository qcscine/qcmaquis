#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "p_dense_matrix/concept/matrix_interface.hpp"
#include "p_dense_matrix/concept/resizable_matrix_interface.hpp"

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

/*BOOST_AUTO_TEST_CASE_TEMPLATE( summ_operation_test, T, test_types )
{
    ambient::layout >> dim(1,1), dim(1,1), dim(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);
    p_dense_matrix<T> d(M_SIZE,M_SIZE);

    c = a + b;
    try{
    c(5,5) = 13.0;
    c(6,5) = 14.0;
    }catch(...){}
    a.remove_rows(128,128);
    b.remove_rows(128,128);
    c.remove_rows(128,128);
    c = a + b;
    a.resize(640,512);
    b.resize(640,512);
    c.resize(640,512);
    c = a + b;
    ambient::playout();
}


BOOST_AUTO_TEST_CASE_TEMPLATE( sql_test, T, test_types )
{
    ambient::layout >> dim(1,1), dim(1,1), dim(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);

    c = a + b;
    c = a + b;
    ambient::playout();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( print_test, T, test_types )
{
    ambient::layout >> dim(1,1), dim(1,1), dim(10,1);

    p_dense_matrix<T> a(8,8);
    ambient::push(ambient::init_double_l,ambient::init_double_c,a); 
    a.remove_rows(0,1);

    std::cout << a;
}
// small matrix class to avoid conflict between p_dense_matrix and dense_matrix //
template<class T>
class Matrix_serial
{
public:
   Matrix_serial<T>(p_dense_matrix<T> &a){
       mx.resize(a.num_rows(), a.num_cols());
       for(int i=0; i<a.num_rows(); i++){
           for(int j=0; j<a.num_cols(); j++){
               mx(i,j) = a(i,j); //very slow .... 
           }
       } 
   } 
  
   Matrix_serial<T>i(Matrix_serial<T> &a){
       mx(a.mx);  
   }

   T& operator()(size_t i, size_t j){ return mx(i,j); }
   boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> mx;
   friend std::ostream& operator<<(std::ostream& os, const Matrix_serial& a){
       for(int i=0; i<a.mx.num_rows(); i++){
           for(int j=0; j<a.mx.num_cols(); j++){
               os<< a(i,j);
               return os; 
           }
       } 
   }    
};




BOOST_AUTO_TEST_CASE_TEMPLATE( SYEV_test, T, test_types ) 
{
    int info; 
    double res = 0; 
    double epsilon = 0.0000000000000001;
 
//first pure ambient run
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    
    p_dense_matrix<T> A(M_SIZE,M_SIZE); 
    p_dense_matrix<T> U(M_SIZE,M_SIZE); 
    p_dense_matrix<T> V(M_SIZE,M_SIZE); 
 
    typename::associated_diagonal_matrix<p_dense_matrix<T> >::type S(M_SIZE,-1); 
 
    A.set_init(ambient::random_i<T>);
    ambient::syev(A,U,S);
    ambient::playout();
 
//second matrix run without ambient

    //  I hope your matrixes are not huges because copy from ambient to serial element by element //
    Matrix_serial<T> A_serial(A);
    Matrix_serial<T> A_serial_bis(A_serial);

    boost::numeric::ublas::vector<T> S_serial_syev(M_SIZE); 
    boost::numeric::ublas::vector<T> S_serial_syevd(M_SIZE); 

    info = boost::numeric::bindings::lapack::syevd('V', A_serial.mx, S_serial_syevd);
    if (info != 0){ throw std::runtime_error("Error in SYEV!");}
    info = boost::numeric::bindings::lapack::syev('V', A_serial_bis.mx, S_serial_syev);
    if (info != 0){ throw std::runtime_error("Error in SYEVD!");}

    for(int i=0 ; i<S_serial_syev.size();i++)
         BOOST_CHECK_CLOSE(S_serial_syev[i],S_serial_syevd[i],epsilon);

    std::reverse(S_serial_syev.begin(), S_serial_syev.end());

    for(int i=0 ; i<S_serial_syev.size();i++){
        res = (fabs(S[i]-S_serial_syev[i]))/fabs(epsilon*S_serial_syev[i]); 
        if(res > 16){ // 16 = Dongarra = number of digit
             printf("validation syev failed, res %.16f Ambient: %.16f Lapack: %.16f \n", res, S_serial_syev[i], S[i]);
        }
     }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SVD_test, T, test_types ) 
{ 
    double res = 0; 
    double epsilon = 0.0000000000000001;
 
//first pure ambient run
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    
    p_dense_matrix<T> A(M_SIZE,N_SIZE); 
    p_dense_matrix<T> U(M_SIZE,M_SIZE); 
    p_dense_matrix<T> V(N_SIZE,N_SIZE); 
 
    typename::associated_diagonal_matrix<p_dense_matrix<T> >::type S(N_SIZE,-1); 
 
    A.set_init(ambient::random_i<T>);
    ambient::svd(A,U,V,S);
    ambient::playout();
 
//second matrix run without ambient

    //  I hope your matrixes are not huges because copy from ambient to serial element by element // 

    Matrix_serial<T> A_serial(A);
    Matrix_serial<T> U_serial(U);
    Matrix_serial<T> V_serial(V);

    boost::numeric::ublas::vector<T> S_serial(std::min(M_SIZE,N_SIZE)); 
 
    int info = boost::numeric::bindings::lapack::gesvd('S', 'S', A_serial.mx, S_serial, U_serial.mx, V_serial.mx);
    if (info != 0){ throw std::runtime_error("Error in SVD!");}
 
    for(int i=0 ; i<S_serial.size();i++){
        res = (fabs(S[i]-S_serial[i]))/fabs(epsilon*S_serial[i]); 
        if(res > 16){ // 16 = Dongarra = number of digit
             printf("validation svd failed, res %.16f Ambient: %.16f Lapack: %.16f \n", res, S[i], S_serial[i]);
        }
     }
}
*/
/*
BOOST_AUTO_TEST_CASE_TEMPLATE( heap_manual_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    p_dense_matrix<T,ambient::MANUAL>* A = new p_dense_matrix<T, ambient::MANUAL>(M_SIZE,M_SIZE);
    p_dense_matrix<T,ambient::MANUAL>* U = new p_dense_matrix<T, ambient::MANUAL>(M_SIZE,M_SIZE); 
    p_dense_matrix<T,ambient::MANUAL>* V = new p_dense_matrix<T, ambient::MANUAL>(M_SIZE,M_SIZE); 
 
    ambient::push(ambient::init_double_l,ambient::init_double_c,*A); 
    ambient::push(ambient::init_double_l,ambient::init_double_c,*U); 
    ambient::push(ambient::init_double_l,ambient::init_double_c,*V); 
    *A = *U + *V;
    ambient::playout();
    delete A; delete U; delete V;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( heap_auto_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    p_dense_matrix<T,ambient::WEAK>* A = new p_dense_matrix<T,ambient::WEAK>(M_SIZE,M_SIZE);
    p_dense_matrix<T,ambient::WEAK>* U = new p_dense_matrix<T,ambient::WEAK>(M_SIZE,M_SIZE); 
    p_dense_matrix<T,ambient::WEAK>* V = new p_dense_matrix<T,ambient::WEAK>(M_SIZE,M_SIZE); 

    ambient::push(ambient::init_double_l,ambient::init_double_c,*A); 
    ambient::push(ambient::init_double_l,ambient::init_double_c,*U); 
    ambient::push(ambient::init_double_l,ambient::init_double_c,*V); 
    *A = *U + *V;
    ambient::playout();
    delete A;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( stack_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    p_dense_matrix<T> A(M_SIZE,M_SIZE);
    p_dense_matrix<T> U(M_SIZE,M_SIZE); 
    p_dense_matrix<T> V(M_SIZE,M_SIZE); 

    ambient::push(ambient::init_double_l,ambient::init_double_c,A); 
    ambient::push(ambient::init_double_l,ambient::init_double_c,U); 
    ambient::push(ambient::init_double_l,ambient::init_double_c,V); 
    A = U + V;
    ambient::playout();
}*/
/*
template<typename T>
void remote_gemm(p_dense_matrix<T> A, p_dense_matrix<T> B, p_dense_matrix<T> C)
{
    C(0,0) = 1;
    C = A * B;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( functionally_remote_gemm_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 
    size_t task_size = 16;

    p_dense_matrix<T> A(task_size,task_size);
    p_dense_matrix<T> B(task_size,task_size); 
    p_dense_matrix<T> C(task_size,task_size); 

    remote_gemm(A, B, C);

    ambient::playout();
    std::cout << C;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( gemm_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 

    size_t task_size = M_SIZE;

    p_dense_matrix<T> A(task_size,task_size);
    p_dense_matrix<T> B(task_size,task_size); 
    p_dense_matrix<T> C(task_size,task_size);

    p_dense_matrix<T> A2(task_size,task_size);
    p_dense_matrix<T> B2(task_size,task_size); 
    p_dense_matrix<T> C2(task_size,task_size); 

    C = A * B;
    C2 = A2 * B2;
    ambient::playout();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 

    p_dense_matrix<T> B(5,6);
    p_dense_matrix<T> BT = blas::transpose(B);

    std::cout << B;
    std::cout << BT;
}

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

BOOST_AUTO_TEST_CASE_TEMPLATE( std_replacement_test, T, test_types ) 
{
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 

    double scut = 0.8;
    size_t nmax = 4;
    size_t* out_value_findif     = (size_t*)malloc(sizeof(size_t));
    double* out_value_accumulate = (double*)malloc(sizeof(double));
    double* out_value_max        = (double*)malloc(sizeof(double));

    typename blas::associated_vector< p_dense_matrix<T> >::type S(M_SIZE,-1);

    S.get_data().set_init(ambient::random_i<T>);

    blas::sort<p_dense_matrix<T> >(S);
    blas::reverse<p_dense_matrix<T> >(S);
    blas::find_if<p_dense_matrix<T> >(S,scut,out_value_findif);
    blas::accumulate<p_dense_matrix<T> >(S,out_value_findif,out_value_accumulate);
    blas::max<p_dense_matrix<T> >(S,scut,nmax,out_value_max); 

    std::cout << S << std::endl;
    printf("out_find_if, out_accumulate, out_max %x %.2f %.2f \n ", *out_value_findif, *out_value_accumulate, *out_value_max); 
    // TODO: CHECK ACCUMULATE
}

BOOST_AUTO_TEST_CASE_TEMPLATE( trace_test, T, test_types ) 
{ 
    ambient::layout >> dim(2,2), dim(2,2), dim(10,1); 

    p_dense_matrix<T> A = p_dense_matrix<T>::identity_matrix(5);
    zout << "Trace: " << blas::trace(A) << "\n";
    std::cout << A;
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


