/*
 *  Const.h
 *  SMART_MATRIX
 *
 *  Created by Tim Ewart on 13.10.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */


extern "C" int dgemm_(char *transa, char *transb, int *m, int *n,
					   int *k, double *alpha, double *a, int *lda,
					   double *b, int *ldb, double *beta, double *c, 
					   int *ldc);


extern "C" void dgesdd_(char *jobz,  int *m,  int  *n, double  * a,  int *lda,double *s,double  *u, int *ldu, 
						 double  *vt,  int  *ldvt, double *work,  int *lwork,  int *iwork, int *info);


extern "C" void dgesvd_(char *jobu, char *jobvt,  int *m,  int *n,
						double *a,  int *lda, double *s, double *u,
						 int *ldu, double *vt,  int *ldvt, double *work,
						 int *lwork,  int *info);


extern "C" void dgeqrf_(int *m,  int *n,
						double *a,  int *lda, double *tau, 
						double *work, int *lwork,  int *info);



extern "C" void dsytrf_(char* uplo, int* n, double *A, int * LDA, int * ipiv, double * work, int * lwork, int * info);



extern "C" int dsyevd_(char *jobz, char *uplo,  int *n,  
					   double *A,  int *lA, double *w, double *work,  int *lwork, 
					    int *iwork,  int *liwork,  int *info);


extern "C" void zgemm_(char *transa, char *transb, int *m, int *n,
					  int *k, dcomplex* alpha, dcomplex*a, int *lda,
					  dcomplex *b, int *ldb, dcomplex *beta, dcomplex *c, 
					  int *ldc);


extern "C" void zgesdd_(char *jobz,  int *m,  int  *n, dcomplex* a,  int *lda,double *s,dcomplex* u, int *ldu, 
						dcomplex  *vt,  int  *ldvt, dcomplex *work,  int *lwork,double *rwork, int *iwork, int *info);


extern "C" void zgesvd_(char *jobu, char *jobvt,  int *m,  int *n,
						dcomplex *a,  int *lda, double *s, dcomplex *u,
						int *ldu, dcomplex *vt,  int *ldvt, dcomplex *work,
						int *lwork, double *rwork,  int *info);


extern "C" void zgeqrf_(int *m,  int *n,
						dcomplex *a,  int *lda, dcomplex *tau, 
						dcomplex *work, int *lwork,  int *info);


extern "C" void cheevd_(char *jobz, char *uplo,  int *n,  
					   fcomplex *A,  int *ldA, float *w, fcomplex *work,  int *lwork, float* rwork, int* lrwork,
					   int *iwork,  int *liwork,  int *info);