
#ifndef M_MATRIX_GPU_FUNCTIONS
#define M_MATRIX_GPU_FUNCTIONS

#include <vector>
#include "dense_matrix/dense_matrix.h"

#include "cuda.h"
//#include "cuda_runtime_api.h" // need for cudaMemcpy, I do not know why !!!!
#include "cublas.h"

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
        
        cublasStatus status;
        
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
        
        int size_matrixA = m*k;
        status = cublasAlloc(size_matrixA, sizeof(double), (void**)&pDataA );
        Check(status , "Set Alloc failed A");
        cublasSetMatrix(m, k, sizeof(double), &lhs(0,0), lda, pDataA, m);	
        Check(status , "Set Matrix failed A");
        
        
        int size_matrixB = k*n_gpu;
        cublasAlloc(size_matrixB, sizeof(double), (void**)&pDataB );
        Check(status , "Set Alloc failed B");
        cublasSetMatrix(k, n_gpu, sizeof(double), &rhs(0,0), ldb, pDataB, k);
        Check(status , "Set Matrix failed B");
        
        int size_matrixC = m*n_gpu;
        cublasAlloc(size_matrixC, sizeof(double), (void**)&pDataC);
        Check(status , "Set Alloc failed C");
        cublasSetMatrix(m, n_gpu, sizeof(double), &res(0,0), ldc, pDataC, m);	
        Check(status , "Set Matrix failed C");
        
        cublasDgemm('n', 'n', m, n_gpu, k, alpha, pDataA, m, pDataB, k, beta, pDataC, m); 
        status = cublasGetError();
        Check(status , "dgemm failed ");
        if (n_cpu > 0)
            dgemm_("n", "n", &m, &n_cpu, &k, &alpha, &lhs(0,0), &lda, &rhs(0,0)+ldb*n_gpu, &ldb, &beta, &res(0,0)+ldc*n_gpu, &ldc); 
        
        status = cublasGetMatrix (m, n_gpu, sizeof(double), pDataC, m, &res(0,0), m); 
        Check(status , "Get data failed");
        
        cublasFree(pDataA);
        cublasFree(pDataB);
        cublasFree(pDataC);
        
    }
    
}
#endif
