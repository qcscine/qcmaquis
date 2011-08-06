#define NODE_COUNT 1

extern "C" {
    void Cblacs_get( int, int, int* );
    void Cblacs_gridinit( int*, const char*, int, int );
    void Cblacs_gridinfo( int, int*, int*, int*, int* );
    void Cblacs_gridexit(int *);
    void Cblacs_exit(int *);
    int  Csys2blacs_handle( MPI_Comm );
    int  numroc_( int*, int*, int*, int*, int* );
    void descinit_( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void pdgemm_(const char*,const char*,int*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
    void pdgesvd_(char *jobu, char *jobvt, int *m, int *n,double *a, int *ia, int *ja, int *desca,double *s,double *u, int *iu, int *ju, int *descu,double *vt, int *ivt, int *jvt, int *descvt,double *work, int *lwork,int *info);
    void pdsyev_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork, int *info);
}

void gemm_c_scalapack(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
#ifdef SCALAPACK
    int info, ictxt, nprow, npcol, myrow, mycol, bn;
    int desca[9], descb[9], descc[9];
    int ZERO=0, ONE=1;
    int anm = get_grid_dim(a).y*get_mem_t_dim(a).y;
    int ann = get_grid_dim(a).x*get_mem_t_dim(a).x; 
    int bnm = get_grid_dim(b).y*get_mem_t_dim(b).y;
    int bnn = get_grid_dim(b).x*get_mem_t_dim(b).x; 
    int cnm = get_grid_dim(c).y*get_mem_t_dim(c).y;
    int cnn = get_grid_dim(c).x*get_mem_t_dim(c).x; 

    nprow = scope.np;
    npcol = scope.nq; 
    bn = get_mem_dim(c).x*get_item_dim(c).x;

    ictxt = Csys2blacs_handle(scope.get_group()->mpi_comm);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int ma = numroc_(&anm, &bn, &myrow, &ZERO, &nprow);
    int na = numroc_(&ann, &bn, &mycol, &ZERO, &npcol);
    int mb = numroc_(&bnm, &bn, &myrow, &ZERO, &nprow);
    int nb = numroc_(&bnn, &bn, &mycol, &ZERO, &npcol);
    int mc = numroc_(&cnm, &bn, &myrow, &ZERO, &nprow);
    int nc = numroc_(&cnn, &bn, &mycol, &ZERO, &npcol);

    descinit_(desca, &anm, &ann, &bn, &bn, &ZERO, &ZERO, &ictxt, &ma, &info);
    descinit_(descb, &bnm, &bnn, &bn, &bn, &ZERO, &ZERO, &ictxt, &mb, &info);
    descinit_(descc, &cnm, &cnn, &bn, &bn, &ZERO, &ZERO, &ictxt, &mc, &info);

    assert(current(a).layout->get_list().size() != 0);
    current(a).solidify(current(a).layout->get_list());
    current(b).solidify(current(b).layout->get_list());
    current(c).solidify(current(c).layout->get_list());

    double alpha = 1.0; double beta = 0.0;
    pdgemm_("N","N",&anm,&bnn,&bnm,&alpha,(double*)current(a).data,&ONE,&ONE,desca,(double*)current(b).data,&ONE,&ONE,descb,&beta,(double*)current(c).data,&ONE,&ONE,descc);
    current(c).disperse(current(c).layout->get_list());
#endif
}

void svd_c_scalapack(const p_dense_matrix<double>& a, int& m, int& n, p_dense_matrix<double>& u, p_dense_matrix<double>& vt, p_dense_matrix<double>& s)
{
#ifndef SCALAPACK
    int info, ictxt, nprow, npcol, myrow, mycol, bn;
    int desca[9], descv[9], descu[9];
    int ZERO=0, ONE=1;
    int k = std::min(n,m);

    nprow = scope.np;
    npcol = scope.nq; 
    bn = get_mem_t_dim(a).x;
    ictxt = Csys2blacs_handle(scope.get_group()->mpi_comm);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int ma = numroc_( &m, &bn, &myrow, &ZERO, &nprow ); // to check // real lda or not?
    int na = numroc_( &n, &bn, &mycol, &ZERO, &npcol );
    int mu = numroc_( &m, &bn, &myrow, &ZERO, &nprow );
    int nu = numroc_( &k, &bn, &mycol, &ZERO, &npcol );
    int mv = numroc_( &k, &bn, &myrow, &ZERO, &nprow );
    int nv = numroc_( &n, &bn, &mycol, &ZERO, &npcol );

    descinit_(desca, &m, &n, &bn, &bn, &ZERO, &ZERO, &ictxt, &ma, &info);
    descinit_(descu, &m, &k, &bn, &bn, &ZERO, &ZERO, &ictxt, &mu, &info);
    descinit_(descv, &k, &n, &bn, &bn, &ZERO, &ZERO, &ictxt, &mv, &info);
   
    assert(current(a).layout->get_list().size() != 0);
    current(a).solidify(current(a).layout->get_list());
    current(u).solidify(current(u).layout->get_list());
    current(vt).solidify(current(vt).layout->get_list());
    current(s).solidify(current(s).layout->get_list());

    int lwork = -1;
    double wkopt;

    //SCALAPACK, first, dry run to allocate buffer
    pdgesvd_("V","V",&m,&n,(double*)breakdown(a).data,&ONE,&ONE,desca,(double*)breakdown(s).data,(double*)breakdown(u).data,&ONE,&ONE,descu,(double*)breakdown(vt).data,&ONE,&ONE,descv,&wkopt,&lwork,&info);

    lwork = static_cast<int>(wkopt);
    double* work = (double*)malloc( lwork*sizeof(double) );

    pdgesvd_("V","V",&m,&n,(double*)breakdown(a).data,&ONE,&ONE,desca,(double*)breakdown(s).data,(double*)breakdown(u).data,&ONE,&ONE,descu,(double*)breakdown(vt).data,&ONE,&ONE,descv,work,&lwork,&info);

    current(u).disperse(current(u).layout->get_list());
    current(vt).disperse(current(vt).layout->get_list());
    current(s).disperse(current(s).layout->get_list());

    free( (void*)work );
#else
/* Locals */
    int lda = get_grid_dim(a).y*get_mem_t_dim(a).y;
    int ldu = get_grid_dim(u).y*get_mem_t_dim(u).y;
    int ldvt = get_grid_dim(vt).y*get_mem_t_dim(vt).y;
    int info, lwork;
    double wkopt;
    double* work;
    assert(current(a).layout->get_list().size() != 0);
    current(a).solidify(current(a).layout->get_list());
    current(u).solidify(current(u).layout->get_list());
    current(vt).solidify(current(vt).layout->get_list());
    current(s).solidify(current(s).layout->get_list());
/* Query and allocate the optimal workspace */
    lwork = -1;
    dgesvd( "S", "S", &m, &n, (double*)breakdown(a).data, &lda, (double*)breakdown(s).data, (double*)breakdown(u).data, &ldu, (double*)breakdown(vt).data, &ldvt, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
/* Compute SVD */
    dgesvd( "S", "S", &m, &n, (double*)breakdown(a).data, &lda, (double*)breakdown(s).data, (double*)breakdown(u).data, &ldu, (double*)breakdown(vt).data, &ldvt, work, &lwork, &info );
/* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm computing SVD failed to converge.\n" );
        exit( 1 );
    }
    current(u).disperse(current(u).layout->get_list());
    current(vt).disperse(current(vt).layout->get_list());
    current(s).disperse(current(s).layout->get_list());
    free( (void*)work );
#endif
}

void syev_c_scalapack(const p_dense_matrix<double>& a, int& m, p_dense_matrix<double>& w)
{
     int lda = get_grid_dim(a).y*get_mem_t_dim(a).y;
     int info, lwork = -1;

     double wkopt;
     double* work;
     assert(current(a).layout->get_list().size() != 0);
     current(a).solidify(current(a).layout->get_list());
     current(w).solidify(current(w).layout->get_list());
     
     dsyev("V","U",&m,(double*)breakdown(a).data,&lda,(double*)breakdown(w).data,&wkopt,&lwork,&info);
     lwork = (int)wkopt;
     work = (double*)malloc( lwork*sizeof(double) );
     dsyev("V","U",&m,(double*)breakdown(a).data,&lda,(double*)breakdown(w).data,work,&lwork,&info);

     if( info > 0 ) {
         printf( "The algorithm computing SYEV failed to converge.\n" );
         exit( 1 );
     }

     current(a).disperse(current(a).layout->get_list());
     current(w).disperse(current(w).layout->get_list());
     free( (void*)work ); 

/*#ifdef SCALAPACK
    int info, ictxt, nprow, npcol, myrow, mycol, bn;
    int desca[9], descz[9];
    int ZERO=0, ONE=1;
    int nm = get_grid_dim(a).y*get_mem_t_dim(a).y;
    int nn = get_grid_dim(a).x*get_mem_t_dim(a).x;

    nprow = scope.np;
    npcol = scope.nq; 
    bn = get_mem_t_dim(a).x;
    ictxt = Csys2blacs_handle(scope.get_group()->mpi_comm);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int ma = numroc_(&nm, &bn, &myrow, &ZERO, &nprow);
    int na = numroc_(&nn, &bn, &mycol, &ZERO, &npcol);
    int mz = numroc_(&nm, &bn, &myrow, &ZERO, &nprow);
    int nz = numroc_(&nn, &bn, &mycol, &ZERO, &npcol);
      cout << ma << " " << na << " " << mz << " " << nz << endl;

    descinit_(desca, &nm, &nn, &bn, &bn, &ZERO, &ZERO, &ictxt, &ma, &info);
    descinit_(descz, &nm, &nn, &bn, &bn, &ZERO, &ZERO, &ictxt, &mz, &info);
   
    assert(current(a).layout->get_list().size() != 0);
    current(a).solidify(current(a).layout->get_list());
    current(w).solidify(current(w).layout->get_list());
    current(z).solidify(current(z).layout->get_list());

    int lwork = -1;
    double wkopt;
      cout << "breakdown " << *(double*)breakdown(a).data << endl;

    pdsyev_("V","U",&nm,(double*)breakdown(a).data,&ONE,&ONE,desca,(double*)breakdown(w).data,(double*)breakdown(z).data,&ONE,&ONE,descz,&wkopt,&lwork,&info);

    lwork = static_cast<int>(wkopt);
    double * work = new double[lwork];

      cout << lwork << endl;
    pdsyev_("V","U",&nm,(double*)breakdown(a).data,&ONE,&ONE,desca,(double*)breakdown(w).data,(double*)breakdown(z).data,&ONE,&ONE,descz,work,&lwork,&info);
    cout << "w " << *(double*)breakdown(w).data << endl;
    cout << "z " << *(double*)breakdown(z).data << endl;
    cout << info << endl;

    current(w).disperse(current(w).layout->get_list());
    current(z).disperse(current(z).layout->get_list());

    delete[] work;
#endif
*/
}

void heev_c_scalapack(const p_dense_matrix<double>& a, int& m, p_dense_matrix<double>& w)
{
     int lda = get_grid_dim(a).y*get_mem_t_dim(a).y;
     int info, lwork = -1;

     double wkopt;
     double* work;
     assert(current(a).layout->get_list().size() != 0);
     current(a).solidify(current(a).layout->get_list());
     current(w).solidify(current(w).layout->get_list());
     
     dsyev("V","U",&m,(double*)breakdown(a).data,&lda,(double*)breakdown(w).data,&wkopt,&lwork,&info);
     lwork = (int)wkopt;
     work = (double*)malloc( lwork*sizeof(double) );
     dsyev("V","U",&m,(double*)breakdown(a).data,&lda,(double*)breakdown(w).data,work,&lwork,&info);

     if( info > 0 ) {
         printf( "The algorithm computing SYEV failed to converge.\n" );
         exit( 1 );
     }

     current(a).disperse(current(a).layout->get_list());
     current(w).disperse(current(w).layout->get_list());
     free( (void*)work );
}
