#include "utils/testing.hpp"

TEST_CASE( "Matrix multiplication performance measured", "[blas::gemm]" )
{
    measurement params;
    size_t x = params.num_cols();
    size_t y = params.num_rows();

    double* ad;
    double* bd;
    double* cd;
    
    int m,n,k;
    int lda,ldb,ldc;
    double alpha(1.0), beta(1.0);

    m = (int)x;
    n = (int)y;
    k = (int)y;

    lda = ldb = ldc = (int)x;

    ad = new double[x*y]; 
    bd = new double[x*y]; 
    cd = new double[x*y]; 
  
    memset((void*)cd,0,x*y*sizeof(double));
   
    for(int i(0); i< x*y ; ++i){
        ad[i] = ambient::utils::Rd();
        bd[i] = ambient::utils::Rd();
    }
       
    measurement::timer time("gemm");
    time.begin();
    ambient::numeric::mkl::blas<double>::gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
    time.end();

    params.report(gflops::gemm, time.get_time());
    delete[] ad;
    delete[] bd;
    delete[] cd;
}
