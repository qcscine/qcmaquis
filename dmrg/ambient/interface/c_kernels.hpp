#define NODE_COUNT 1

#include "mkl.h"

extern "C" {
    void Cblacs_get( int, int, int* );
    void Cblacs_gridinit( int*, const char*, int, int );
    void Cblacs_gridinfo( int, int*, int*, int*, int* );
    int Csys2blacs_handle( MPI_Comm );
    int numroc_( int*, int*, int*, int*, int* );
    void descinit_( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void pdgemm_(const char*,const char*,int*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
}

/*
 --- --- ---       --- --- ---       --- --- ---
| 0 |   |   |     | 0 | ! | ! |     | 0 |   |   |
 --- --- ---       --- --- ---       --- --- ---
| 0 |   |   |  x  |   |   |   |  =  | 0 |   |   |
 --- --- ---       --- --- ---       --- --- ---
| 0 |   |   |     |   |   |   |     | 0 |   |   |
 --- --- ---       --- --- ---       --- --- ---
partial reduce?
*/

void gemm_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c){
// todo
    double* ad = breakdown(a)(get_group_id(c).x, get_group_id(c).y);

// ok, let's try to use MKL
    int m = 3;
    int n = 3;
    int k = 3;
    double* a_ = (double*)malloc(sizeof(double)*m*n);
    int lda = 3;
    double* b_ = (double*)malloc(sizeof(double)*m*n);
    int ldb = 3;
    double* c_ = (double*)malloc(sizeof(double)*m*n);
    int ldc = 3;
    double alpha = 1;
    double beta = 0;

    for(int j=0; j<n; j++)
    for(int i=0; i<m; i++){
        a_[j*lda+i] = i+j;
        b_[j*lda+i] = i+j;
        c_[j*lda+i] = i+j;
    }
    if(scope.get_rank() == 0){
      for(int i=0; i<m; i++){
        for(int j=0; j<n; j++) printf("%.2f	", a_[j*lda+i]);
        printf("\n");
      }
      printf("\n");
    }

    dgemm("N","N", &m, &n, &k, &alpha, a_, &lda, b_, &ldb, &beta, c_, &ldc);

    if(scope.get_rank() == 0){
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++) printf("%.2f	", c_[j*ldc+i]);
        printf("\n");
    }
    printf("\n\nDone!\n\n"); }
}

void pdgemm_c_kernel(p_dense_matrix<double>& a, p_dense_matrix<double>& b, p_dense_matrix<double>& c){
//    printf("R%d: Executing ScaLAPACK PDGEMM kernel\n", scope.get_rank());

/*    int i, j, k;
    int bhandle, ictxt, nprow, npcol, myrow, mycol,nb;
    nprow = NODE_COUNT; npcol = (int)(scope.get_size()/NODE_COUNT); 
    nb = breakdown(c).get_group_dim().x*breakdown(c).get_item_dim().x;
    int M=breakdown(c).get_grid_dim()*breakdown(c).get_group_dim().x*breakdown(c).get_item_dim().x;
    int info,itemp;
    int ZERO=0,ONE=1;
 
    bhandle = Csys2blacs_handle(scope.get_group()->mpi_comm);
    ictxt = bhandle;
    //Cblacs_get( -1, 0, &ictxt );
    Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
    //printf("Number of rows: %d\nNumber of cols: %d\nMy row: %d\nMy col: %d\n", nprow, npcol, myrow, mycol);

    int descA[9],descB[9],descC[9];
    int mA = numroc_( &M, &nb, &myrow, &ZERO, &nprow );
    int nA = numroc_( &M, &nb, &mycol, &ZERO, &npcol );
    int mB = numroc_( &M, &nb, &myrow, &ZERO, &nprow );
    int nB = numroc_( &M, &nb, &mycol, &ZERO, &npcol );
    int mC = numroc_( &M, &nb, &myrow, &ZERO, &nprow );
    int nC = numroc_( &M, &nb, &mycol, &ZERO, &npcol );
    descinit_(descA, &M,   &M,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);
    descinit_(descB, &M,   &M,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mB,  &info);
    descinit_(descC, &M,   &M,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mC,  &info);
    
   // breakdown(a).solidify();
   // breakdown(b).solidify();
   // breakdown(c).solidify();
    double *A = (double*) malloc(mA*nA*sizeof(double));
    double *B = (double*) malloc(mB*nB*sizeof(double));
    double *C = (double*) malloc(mC*nC*sizeof(double));
    for(i=0;i<mA;i++) for(j=0;j<nA;j++){
                 A[j*mA+i]=0.01;
         }
    for(i=0;i<mB;i++) for(j=0;j<nB;j++){
                 B[j*mB+i]=0.01;
         }
    for(i=0;i<mC;i++) for(j=0;j<nC;j++){
                 C[j*mC+i]=0.01;
         }
    printf("M: %d; mA: %d; nA: %d\n", M, mA, nA);
    double alpha = 1.0; double beta = 0.0;
    
//    for(i=0;i<VECTOR_SIZE;i++)
    //pdgemm_("N","N",&M,&M,&M,&alpha,(double*)breakdown(a).data,&ONE,&ONE,descA,(double*)breakdown(b).data,&ONE,&ONE,descB,&beta,(double*)breakdown(c).data,&ONE,&ONE,descC);
    pdgemm_("N","N",&M,&M,&M,&alpha,A,&ONE,&ONE,descA,B,&ONE,&ONE,descB,&beta,C,&ONE,&ONE,descC);
*/
}

void null_c_kernel(p_dense_matrix<double>& a){}
void add_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c){}
void sub_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c){}
void scale_c_kernel(const p_dense_matrix<double>& a, const double& alfa, pinned p_dense_matrix<double>& out){}

/////////////////////
// testing kernels // 

void single_integer_c_kernel(int& input){
    input += 13;
    zout << "single integer kernel: output is " << input << "\n";
}
