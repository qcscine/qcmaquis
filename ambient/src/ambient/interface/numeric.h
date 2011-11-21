/*
* Fortran declaration : native compatibility all math libs (MKL, Goto, ACML, ESSLSMP, ...)
*/

extern "C" {
    double sqrt(double);
    double fabs(double);

    void dgemm_(const char *transa, const char *transb, const int  *m, const int *n, const int *k,
                const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,                                                   
                const double *beta, double *c, const int *ldc);

    void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);

    void dgesvd_( const char* jobu, const char* jobvt, const int* m, 
                  const int* n, double* a, const int* lda, double* s, 
                  double* u, const int* ldu, double* vt, const int* ldvt, 
                  double* work, const int* lwork, int* info );

    void dsyev_( const char* jobz, const char* uplo, const int* n, double* a, 
                 const int* lda, double* w, double* work, const int* lwork, 
                 int* info );
}
