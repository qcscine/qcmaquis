#include "utils/testing.hpp"

#define AA(i,j) AA[(i)*N+(j)]
#define BB(i,j) BB[(i)*N+(j)]

extern "C" void Cblacs_pinfo( int* mypnum, int* nprocs); 
extern "C" void Cblacs_gridinfo( int nContxt, int* nRows, int* nCols, int* nMyRow, int* nMyCol);
extern "C" void Cblacs_gridinit(  int* nContxt, char * order,  int np_row,  int np_col);
extern "C" void Cblacs_gridexit( int nContxt);
extern "C" void Cblacs_exit( int nContxt);
extern "C" void Cblacs_get(int Contxt, int what, int* val); 
extern "C" int numroc_(int* nN, int* nB, int* nIProc, int* nISrcProc, int* nNProcs);
extern "C" void descinit_(int *, int*, int* , int *, int *, int *, int *, int *, int *, int*);
extern "C" void pdgemm_(char *jobu, char *jobvt,
		int *m, int *n, int *k,
		double *alpha, double * a, int *ia, int *ja, int *desca,
		               double * b, int *ib, int *jb, int *descb,
		double *beta , double * c, int *ic, int *jc, int *descc);


TEST_CASE( "Matrix multiplication performance measured", "[pblas::gemm]" )
{
   measurement params;
   int N = params.num_cols();
   int argc=0;
   char ** argv;
   srand(3);
   int i, j, k;
   int minusone = -1;
   int myrank_mpi, nprocs_mpi;
   int ictxt, nprow, npcol, myrow, mycol,nb;

   int info,itemp;
   int zero = 0, one = 1;

   double* AA = (double*)malloc(N*N*sizeof(double));
   double* BB = (double*)malloc(N*N*sizeof(double));
   double* CC = (double*)malloc(N*N*sizeof(double));

   int descA[9],descB[9],descC[9];
   for(i = 0; i < N; i++)
     for(j = 0; j < N; j++){
        AA[i*N + j] = rand();
        BB[i*N + j] = rand();
        CC[i*N + j] = -1;
     }
   measurement::timer time("gemm");
   time.begin();
   Cblacs_pinfo( &myrank_mpi, &nprocs_mpi) ;

   nprow = 1; // 1-d = one row
   npcol = nprocs_mpi; // number of col for 1d cyclic distribution
   nb = 128; // # of element (width x direction of each column) also workgroup into ambient

   Cblacs_get(minusone, zero, &ictxt);
   Cblacs_gridinit( &ictxt, "C", nprow, npcol );
   Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

   int mA = numroc_( &N, &nb, &mycol, &zero, &npcol);
   int mB = numroc_( &N, &nb, &mycol, &zero, &npcol);
   int mC = numroc_( &N, &nb, &mycol, &zero, &npcol);

   //arg 3 = N  number of row of the block issue from the block cyclic decomposition
   //arg 4 = nb number of col of the block issue from the block cyclic decomposition
   //arg 8 = lda issue from numroc here use less because 1d distribution 
 
   descinit_(descA, &N, &N, &nb, &nb, &zero, &zero, &ictxt, &N, &info);
   descinit_(descB, &N, &N, &nb, &nb, &zero, &zero, &ictxt, &N, &info);
   descinit_(descC, &N, &N, &nb, &nb, &zero, &zero, &ictxt, &N, &info);

   double alpha = 1.0; double beta = 0.0;
   pdgemm_("N","N",&N,&N,&N,&alpha,AA,&one,&one,descA,BB,&one,&one,descB,&beta,CC,&one,&one,descC);
   time.end();

   if(myrank_mpi == 0) params.report(gflops::gemm, time.get_time());
   MPI_Barrier(MPI_COMM_WORLD);
   Cblacs_exit(0);
}
