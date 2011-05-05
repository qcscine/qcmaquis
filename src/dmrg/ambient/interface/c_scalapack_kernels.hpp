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

void gemm_c_scalapack_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
#ifdef SCALAPACK
    int info, ictxt, nprow, npcol, myrow, mycol, bn;
    int desca[9], descb[9], descc[9];
    int ZERO=0, ONE=1;
    int nn = get_grid_dim(c).x*get_mem_dim(c).x*get_item_dim(c).x; 
    int nm = get_grid_dim(c).y*get_mem_dim(c).y*get_item_dim(c).y;

    nprow = scope.np;
    npcol = scope.nq; 
    bn = get_mem_dim(c).x*get_item_dim(c).x;

    ictxt = Csys2blacs_handle(scope.get_group()->mpi_comm);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int ma = numroc_(&nm, &bn, &myrow, &ZERO, &nprow);
    int na = numroc_(&nn, &bn, &mycol, &ZERO, &npcol);
    int mb = numroc_(&nm, &bn, &myrow, &ZERO, &nprow);
    int nb = numroc_(&nn, &bn, &mycol, &ZERO, &npcol);
    int mc = numroc_(&nm, &bn, &myrow, &ZERO, &nprow);
    int nc = numroc_(&nn, &bn, &mycol, &ZERO, &npcol);

    descinit_(desca, &nm, &nn, &bn, &bn, &ZERO, &ZERO, &ictxt, &ma, &info);
    descinit_(descb, &nm, &nn, &bn, &bn, &ZERO, &ZERO, &ictxt, &mb, &info);
    descinit_(descc, &nm, &nn, &bn, &bn, &ZERO, &ZERO, &ictxt, &mc, &info);

    assert(current(a).layout->get_list().size() != 0);
    current(a).solidify(current(a).layout->get_list());
    current(b).solidify(current(b).layout->get_list());
    current(c).solidify(current(c).layout->get_list());

    double alpha = 1.0; double beta = 0.0;
    pdgemm_("N","N",&nm,&nn,&nm,&alpha,(double*)current(a).data,&ONE,&ONE,desca,(double*)current(b).data,&ONE,&ONE,descb,&beta,(double*)current(c).data,&ONE,&ONE,descc);
    current(c).disperse(current(c).layout->get_list());
#endif
}

void svd_c_scalapack_kernel(const p_dense_matrix<double>& a, p_dense_matrix<double>& u, p_dense_matrix<double>& v, p_dense_matrix<double>& s)
{
#ifdef SCALAPACK
    int info, ictxt, nprow, npcol, myrow, mycol, bn;
    int desca[9], descv[9], descu[9];
    int ZERO=0, ONE=1;
    int n = get_grid_dim(a).x*get_mem_t_dim(a).x; 
    int m = get_grid_dim(a).y*get_mem_t_dim(a).y;
    int k = std::min(n,m);

    nprow = scope.np;
    npcol = scope.nq; 
    bn = get_mem_t_dim(a).x;
    ictxt = Csys2blacs_handle(scope.get_group()->mpi_comm);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int ma = numroc_( &m, &bn, &myrow, &ZERO, &nprow ); // to check
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
    current(v).solidify(current(v).layout->get_list());
    current(s).solidify(current(s).layout->get_list());

    int lwork = -1;
    double wkopt;

    //SCALAPACK, first, dry run to allocate buffer
    pdgesvd_("V","V",&m,&n,(double*)breakdown(a).data,&ONE,&ONE,desca,(double*)breakdown(s).data,(double*)breakdown(u).data,&ONE,&ONE,descu,(double*)breakdown(v).data,&ONE,&ONE,descv,&wkopt,&lwork,&info);

    lwork = static_cast<int>(wkopt);
    double* work = (double*)malloc( lwork*sizeof(double) );

    pdgesvd_("V","V",&m,&n,(double*)breakdown(a).data,&ONE,&ONE,desca,(double*)breakdown(s).data,(double*)breakdown(u).data,&ONE,&ONE,descu,(double*)breakdown(v).data,&ONE,&ONE,descv,work,&lwork,&info);

    current(u).disperse(current(u).layout->get_list());
    current(v).disperse(current(v).layout->get_list());
    current(s).disperse(current(s).layout->get_list());

    free( (void*)work );
#else
/* Locals */
    int m = get_grid_dim(a).y*get_mem_t_dim(a).y;
    int n = get_grid_dim(a).x*get_mem_t_dim(a).x;
    int lda = m, ldu = m, ldvt = std::min(m,n), info, lwork;
    double wkopt;
    double* work;
    assert(current(a).layout->get_list().size() != 0);
    current(a).solidify(current(a).layout->get_list());
    current(u).solidify(current(u).layout->get_list());
    current(v).solidify(current(v).layout->get_list());
    current(s).solidify(current(s).layout->get_list());
/* Query and allocate the optimal workspace */
    lwork = -1;
    dgesvd( "S", "S", &m, &n, (double*)breakdown(a).data, &lda, (double*)breakdown(s).data, (double*)breakdown(u).data, &ldu, (double*)breakdown(v).data, &ldvt, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
/* Compute SVD */
    dgesvd( "S", "S", &m, &n, (double*)breakdown(a).data, &lda, (double*)breakdown(s).data, (double*)breakdown(u).data, &ldu, (double*)breakdown(v).data, &ldvt, work, &lwork, &info );
/* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm computing SVD failed to converge.\n" );
        exit( 1 );
    }
    current(u).disperse(current(u).layout->get_list());
    current(v).disperse(current(v).layout->get_list());
    current(s).disperse(current(s).layout->get_list());
    free( (void*)work );
#endif
}


void syev_c_scalapack_kernel(const p_dense_matrix<double>& a, p_dense_matrix<double>& w, p_dense_matrix<double>& z)
{
#ifdef SCALAPACK
    int info, ictxt, nprow, npcol, myrow, mycol, bn;
    int desca[9], descw[9], descz[9];
    int ZERO=0, ONE=1;
    int nm = get_grid_dim(a).y*get_mem_t_dim(a).y;
    int nn = get_grid_dim(a).x*get_mem_t_dim(a).x;

    nprow = scope.np;
    npcol = scope.nq; 
    bn = get_mem_dim(a).x*get_item_dim(a).x;
    ictxt = Csys2blacs_handle(scope.get_group()->mpi_comm);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int ma = numroc_(&nm, &bn, &myrow, &ZERO, &nprow);
    int na = numroc_(&nn, &bn, &mycol, &ZERO, &npcol);
    int mz = numroc_(&nm, &bn, &myrow, &ZERO, &nprow);
    int nz = numroc_(&nn, &bn, &mycol, &ZERO, &npcol);
    int mw = numroc_(&nm, &bn, &myrow, &ZERO, &nprow);
    int nw = numroc_(&nn, &bn, &mycol, &ZERO, &npcol);

    descinit_(desca, &nm, &nn, &bn, &bn, &ZERO, &ZERO, &ictxt, &ma, &info);
    descinit_(descz, &nm, &nn, &bn, &bn, &ZERO, &ZERO, &ictxt, &mz, &info);
    descinit_(descw, &nm, &nn, &bn, &bn, &ZERO, &ZERO, &ictxt, &mw, &info);
   
    assert(current(a).layout->get_list().size() != 0);
    current(a).solidify(current(a).layout->get_list());
    current(w).solidify(current(w).layout->get_list());
    current(z).solidify(current(z).layout->get_list());

    int lwork = -1;
    double wkopt;

    /* SCALAPACK, first, dry run to allocate buffer */
    pdsyev_("V","U",&nm,(double*)breakdown(a).data,&ONE,&ONE,desca,(double*)breakdown(z).data,(double*)breakdown(w).data,&ONE,&ONE,descz,&wkopt,&lwork,&info);

    lwork = static_cast<int>(wkopt);
    double* work = (double*)malloc(sizeof(double)*lwork);

    pdsyev_("V","U",&nm,(double*)breakdown(a).data,&ONE,&ONE,desca,(double*)breakdown(z).data,(double*)breakdown(w).data,&ONE,&ONE,descz,work,&lwork,&info);

    current(w).disperse(current(w).layout->get_list());
    current(z).disperse(current(z).layout->get_list());

    free((void*)work); /* clean the working buffer */
#endif
}

