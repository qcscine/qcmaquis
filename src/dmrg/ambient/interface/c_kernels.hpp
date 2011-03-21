#define NODE_COUNT 1

void add_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
    void_pt& profile = breakdown(out);
    double* ad = breakdown(a)(get_group_id(out).x, get_group_id(out).y);
    double* bd = breakdown(b)(get_group_id(out).x, get_group_id(out).y);
    int size = get_group_dim(out).x*get_item_dim(out).x*
               get_group_dim(out).y*get_item_dim(out).y;
//    printf("R%d: Executing plus computation kernel (%d ops)... for out grp %d %d\n", scope.get_rank(), size, get_group_id(out).x, get_group_id(out).y);
//    for(int i=0; i < size; i++){
//        output[i] = ad[i]+bd[i];
//    }

}

void sub_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
// todo
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

void gemm_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& out){
// todo
    double a = breakdown(a).(get_group_id(out).x, get_group_id(out).y);
    printf("I was called actually by %d\n", ambient::rank());

}

void scale_c_kernel(const p_dense_matrix<double>& a, const double& b, pinned p_dense_matrix<double>& out){
// todo
}

void null_c_kernel(p_dense_matrix<double>& a){
    printf("R%d: Executing NULL kernel\n", scope.get_rank());
}


extern "C" {
    void Cblacs_get( int, int, int* );
    void Cblacs_gridinit( int*, const char*, int, int );
    void Cblacs_gridinfo( int, int*, int*, int*, int* );
    int Csys2blacs_handle( MPI_Comm );
    int numroc_( int*, int*, int*, int*, int* );
    void descinit_( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void pdgemm_(const char*,const char*,int*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
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

void single_integer_c_kernel(int& input){
    input += 13;
    zout << "single integer kernel: output is " << input << "\n";
}
