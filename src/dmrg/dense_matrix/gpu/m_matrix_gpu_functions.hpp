
#ifndef M_MATRIX_GPU_FUNCTIONS
#define M_MATRIX_GPU_FUNCTIONS

#include <vector>
#include "dense_matrix/dense_matrix.h"

#include "cuda.h"
//#include "cuda_runtime_api.h" // needed for streams
#include "cublas.h"

#include <utils/timings.h>

namespace gpu {
    
    void Check(cublasStatus const  status, char const * Message)
    {
        if (status != CUBLAS_STATUS_SUCCESS) 
        {
            printf ("%s \n", Message);
            //		cublasShutdown(); 
        }	
        
        if (status == CUBLAS_STATUS_NOT_INITIALIZED) 
        {
            printf ("cublas not init \n");
        }	
        
        if (status == CUBLAS_STATUS_ALLOC_FAILED) 
        {
            printf ("cublas  alloc failed \n ");
        }	

        if (status == CUBLAS_STATUS_MAPPING_ERROR) 
        {
            printf ("cublas  mapping error \n ");
        }	
        
        if (status == CUBLAS_STATUS_INVALID_VALUE) 
        {
            printf ("cublas invalid value \n ");
        }	
        
        if (status == CUBLAS_STATUS_EXECUTION_FAILED) 
        {
            printf ("cublas execution failed \n ");
        }	
        
        
    }
    
    
    void matrix_matrix_multiply(blas::dense_matrix<double,std::vector<double, std::allocator<double> > > const & lhs,
                                blas::dense_matrix<double,std::vector<double, std::allocator<double> > > const & rhs,
                                blas::dense_matrix<double,std::vector<double, std::allocator<double> > >  & res,
                                double ratio=0.)
    {
        /** ration value
         0.0 = full gpu
         1 = ful cpu
         fermi = 0.2
         tesla = 0.4
         */
        
        assert( lhs.num_cols() == rhs.num_rows() );
        
        static Timer timer("gpu gemm");
        timer.begin();
        
        cublasStatus status;
//        cudaStream_t stream[3];
//        for (int i = 0; i < 3; ++i)
//            cudaStreamCreate(&stream[i]);
        
        double * pDataA = 0;
        double * pDataB = 0;
        double * pDataC = 0;
		        
        int n = rhs.num_cols(); 
        int n_cpu = static_cast<int>(n * ratio);
        int n_gpu = n - n_cpu ;
        int m = lhs.num_rows(); 
        int k = lhs.num_cols();
        double alpha = 1.;
        double beta = 0.;
        res.resize(m, n);
        
        int lda = lhs.stride2();
        int ldb = rhs.stride2();
        int ldc = res.stride2();

#ifdef TIME_GEMM_GPU_INTERNAL
        static Timer tA("gpu gemm: A transfer");
        tA.begin();
#endif
        int size_matrixA = m*k;
        status = cublasAlloc(size_matrixA, sizeof(double), (void**)&pDataA );
        Check(status , "Set Alloc failed A");
        status = cublasSetMatrix(m, k, sizeof(double), &lhs(0,0), lda, pDataA, m);	
//        status = cublasSetMatrixAsync(m, k, sizeof(double), &lhs(0,0), lda, pDataA, m, stream[0]);	
        Check(status , "Set Matrix failed A");
        
#ifdef TIME_GEMM_GPU_INTERNAL
        tA.end();
        static Timer tB("gpu gemm: B transfer");
        tB.begin();
#endif
        int size_matrixB = k*n_gpu;
        status = cublasAlloc(size_matrixB, sizeof(double), (void**)&pDataB );
        Check(status , "Set Alloc failed B");
//        status = cublasSetMatrixAsync(k, n_gpu, sizeof(double), &rhs(0,0), ldb, pDataB, k, stream[1]);
        status = cublasSetMatrix(k, n_gpu, sizeof(double), &rhs(0,0), ldb, pDataB, k);
        Check(status , "Set Matrix failed B");
        
#ifdef TIME_GEMM_GPU_INTERNAL
        tB.end();
        static Timer tC("gpu gemm: C allocation");
        tC.begin();
#endif
        int size_matrixC = m*n_gpu;
        status = cublasAlloc(size_matrixC, sizeof(double), (void**)&pDataC);
        Check(status , "Set Alloc failed C");
//        status = cublasSetMatrix(m, n_gpu, sizeof(double), &res(0,0), ldc, pDataC, m);	
//        status = cublasSetMatrixAsync(m, n_gpu, sizeof(double), &res(0,0), ldc, pDataC, m, stream[2]);	
//        Check(status , "Set Matrix failed C");
        
#ifdef TIME_GEMM_GPU_INTERNAL
        tC.end();
        static Timer t1("gpu gemm: gemm call");
        t1.begin();
#endif
        
//        for (int i = 0; i < 3; ++i) {
//            cudaStreamDestroy(stream[i]);
//        }
        
        cublasDgemm('n', 'n', m, n_gpu, k, alpha, pDataA, m, pDataB, k, beta, pDataC, m); 
        status = cublasGetError();
        Check(status , "dgemm failed ");
        if (n_cpu > 0)
            dgemm_("n", "n", &m, &n_cpu, &k, &alpha, &lhs(0,0), &lda, &rhs(0,0)+ldb*n_gpu, &ldb, &beta, &res(0,0)+ldc*n_gpu, &ldc); 
        
        status = cublasGetMatrix (m, n_gpu, sizeof(double), pDataC, m, &res(0,0), ldc); 
        Check(status , "Get data failed");
        
#ifdef TIME_GEMM_GPU_INTERNAL
        t1.end();
#endif
        
        
        cublasFree(pDataA);
        cublasFree(pDataB);
        cublasFree(pDataC);
        
        timer.end();
        
    }
    
}
#endif
