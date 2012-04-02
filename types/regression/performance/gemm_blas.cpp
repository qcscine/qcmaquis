#define BOOST_TEST_MODULE benchmarks
#include <mpi.h>
#include <omp.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>
#include <cmath>

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/algorithms.hpp"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"

#include "types/utils/matrix_cast.h"
#include "utils/timings.h"

#include "utilities.h"

#include "ambient/utils/numeric.h" 

BOOST_AUTO_TEST_CASE_TEMPLATE(gemm_blas, T, test_types)
{
    omp_set_num_threads(T::ValueThread);
    typename T::value_type*ad;
    typename T::value_type*bd;
    typename T::value_type*cd;
    
    int m,n,k;
    int lda,ldb,ldc;
    double alpha(1.0), beta(1.0);

    m =(int) T::ValueX;
    n =(int) T::ValueY;
    k =(int) T::ValueY;

    lda = ldb = ldc = T::ValueX;

    ad = new typename T::value_type[T::ValueX*T::ValueY]; 
    bd = new typename T::value_type[T::ValueX*T::ValueY]; 
    cd = new typename T::value_type[T::ValueX*T::ValueY]; 
  
    memset((void*)cd,0,T::ValueX*T::ValueY*sizeof(typename T::value_type));
   
    for(int i(0); i< T::ValueX*T::ValueY ; ++i){
        ad[i] = Rd();
        bd[i] = Rd();
    }
       
    TimerPTH t1("blas");
    t1.begin();
    gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
    t1.end();

    double time = t1.GetTime();
    double gfl  = GFlopsGemm(n,m,k,time);
    save("TimeGemmBlas.txt",t1,gfl,T::ValueX,T::ValueY,T::ValueThread); 
  
    delete [] ad;   
    delete [] bd;   
    delete [] cd;   
}

