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
| 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
 --- --- ---       --- --- ---       --- --- ---
| 0 | 1 | 2 |  x  | 0 | 1 | 2 |  =  | 0 | 1 | 2 |
 --- --- ---       --- --- ---       --- --- ---
| 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
 --- --- ---       --- --- ---       --- --- ---
partial reduce?
*/

void gemm_c_kernel(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
    int m   = get_group_dim(a).y*get_item_dim(a).y;
    int n   = get_group_dim(b).x*get_item_dim(b).x;
    int k   = get_group_dim(b).y*get_item_dim(b).y;
    int lda = m;
    int ldb = k;
    int ldc = m;
    double alpha = 1.0; 
    double beta  = 0.0;
// a(i,j) => b(j,i) x a(z,j) where z : [0,m)
// current group of matrix a:
    int i = get_group_id(a).y;
    int j = get_group_id(a).x;
// taking (j,i) of b:
    double* bd = current(b)(j,i); // remote
// multiplying with column of a:
    for(int z = 0; z < get_grid_dim(a).y; z++){
        double* ad = current(a)(z,j);
        double* cd = reduced<'+'>(c)(z,i); // a(z,j) x b(j,i) => c(z,i)
        //printf("Performing DGEMM for %d,%d and %d,%d\n", z,j,j,i);
        dgemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
        //for(int ii=0; ii < m; ii++){
        //  for(int jj=0; jj < n; jj++)
        //  printf("%.2f	", cd[jj*ldc + ii]);
        //  printf("\n");
        //}
    }
}


void pdgemm_c_kernel(p_dense_matrix<double>& a, p_dense_matrix<double>& b, p_dense_matrix<double>& c){
//    printf("R%d: Executing ScaLAPACK PDGEMM kernel\n", scope.get_rank());

/*    int i, j, k;
    int bhandle, ictxt, nprow, npcol, myrow, mycol,nb;
    nprow = NODE_COUNT; npcol = (int)(scope.get_size()/NODE_COUNT); 
    nb = get_group_dim(c).x*get_item_dim(c).x;
    int M = get_grid_dim(c)*get_group_dim(c).x*get_item_dim(c).x;
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


void copy_c_kernel(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a)
{    
    int i = get_group_id(a).y;
    int j = get_group_id(a).x;
    double* a_elements  = current(a)(i,j);
    double* ac_elements = current(ac)(i,j);
    memcpy(ac_elements, a_elements, sizeof(double)*(get_group_dim(a)*get_item_dim(a)));
}

void remove_rows_c_kernel(pinned p_dense_matrix<double>& a, const size_t& i_mark, const size_t& k)
{
    size_t i   = get_group_id(a).y;
    size_t j   = get_group_id(a).x;
    size_t lda = get_group_t_dim(a).y;

    size_t remains_u = i_mark % lda;
    size_t remains_l = lda - (remains_u+k) % lda;
    size_t shift     = k / lda + (k % lda ? 1 : 0);

    if(remains_u + remains_l < lda){
// get two following blocks
        //get: i+shift-1 
        //get: i+shift

    }else{
// get only one following block
        //get: i+shift;

    }
/*
    size_t group_i_mark         = i_mark / lda;
    size_t group_element_i_mark = i_mark % lda;
    size_t group_rows_to_replace = lda - (starting_delete_element_i+1);

    size_t ending_delete_group_i   = (i_mark+k) / lda;
    size_t ending_delete_element_i = (i_mark+k) % lda;
    size_t group_rows_to_replace_w = lda - ending_delete_element_i;

    int movement_lda = min(k, (lda-starting_delete_element_i));

    if(i == starting_delete_group_i){
        // proceed with shifting
        replacement = current(a)(ending_delete_group_i,j);
        if(group_rows_to_replace > group_rows_to_replace_w)
        replacement = current(a)(ending_delete_group_i + 1,j);
    }else if(i > group_start_i){
        // proceed with shifting
        replacement = current(a)(ending_delete_group_i + (i - starting_delete_group_i),j);
        if(group_rows_to_replace > group_rows_to_replace_w)
        replacement = current(a)(ending_delete_group_i + (i - starting_delete_group_i) + 1,j);
        // move my data first 
    }

    int element_i = i % (this->profile->get_group_dim().y*this->profile->get_item_dim().y);
       // for(size_type j = 0; j < this->cols; ++j) // for each column, copy the rows > i+k to k rows  up
       //     memmove(&this->data[this->lda*j + i], &this->data[this->lda*j + i + k], sizeof(T)*(this->rows-i-k));
*/
}

void remove_cols_c_kernel(pinned p_dense_matrix<double>& a, const size_t& j_mark, const size_t& k)
{

}

void null_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c){ }
void add_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c){ printf("Executed add kernel\n"); }
void sub_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c){}
void scale_c_kernel(const p_dense_matrix<double>& a, const double& alfa, pinned p_dense_matrix<double>& out){}

/////////////////////
// testing kernels // 

void single_integer_c_kernel(int& input){
    input += 13;
    zout << "single integer kernel: output is " << input << "\n";
}
