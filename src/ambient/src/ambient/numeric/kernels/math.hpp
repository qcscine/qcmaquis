#ifndef AMBIENT_NUMERIC_MATRIX_KERNELS_MATH
#define AMBIENT_NUMERIC_MATRIX_KERNELS_MATH

#define PLASMA_IB        64
#define PlasmaNoTrans    111
#define PlasmaTrans      112
#define PlasmaUpper      121
#define PlasmaLower      122
#define PlasmaUpperLower 123
#define PlasmaLeft       141
#define PlasmaRight      142
#define PlasmaForward    391
#define PlasmaColumnwise 401

typedef int PLASMA_enum;

extern "C" {

    int LAPACKE_dbdsdc(int matrix_order, char uplo, char compq, int n, double* d, double* e, double* u, int ldu, double* vt, int ldvt, double* q, int* iq);
    double     ddot_(const int*, const double*, const int*, 
                     const double*, const int*);
    void     dgebd2_(int*, int*, double*, 
                     int*, double*, double*, 
                     double*, double*, double*, int*);
    int  CORE_dgeqrt(int M, int N, int IB,
                     double *A, int LDA,
                     double *T, int LDT,
                     double *TAU, double *WORK);
    int  CORE_dormqr(PLASMA_enum side, PLASMA_enum trans,
                     int M, int N, int K, int IB,
                     const double *V, int LDV,
                     const double *T, int LDT,
                     double *C, int LDC,
                     double *WORK, int LDWORK);
    int  CORE_dtsqrt(int M, int N, int IB,
                     double *A1, int LDA1,
                     double *A2, int LDA2,
                     double *T, int LDT,
                     double *TAU, double *WORK);
    int  CORE_dtsmqr(PLASMA_enum side, PLASMA_enum trans,
                     int M1, int N1, int M2, int N2, int K, int IB,
                     double *A1, int LDA1,
                     double *A2, int LDA2,
                     const double *V, int LDV,
                     const double *T, int LDT,
                     double *WORK, int LDWORK);
    int  CORE_dgelqt(int M, int N, int IB,
                     double *A, int LDA,
                     double *T, int LDT,
                     double *TAU,
                     double *WORK);
    int  CORE_dormlq(PLASMA_enum side, PLASMA_enum trans,
                     int M, int N, int IB, int K,
                     const double *V, int LDV,
                     const double *T, int LDT,
                     double *C, int LDC,
                     double *WORK, int LDWORK);
    int  CORE_dtslqt(int M, int N, int IB,
                     double *A1, int LDA1,
                     double *A2, int LDA2,
                     double *T, int LDT,
                     double *TAU, double *WORK);
    int  CORE_dtsmlq(PLASMA_enum side, PLASMA_enum trans,
                     int M1, int N1, int M2, int N2, int K, int IB,
                     double *A1, int LDA1,
                     double *A2, int LDA2,
                     const double *V, int LDV,
                     const double *T, int LDT,
                     double *WORK, int LDWORK);
    void CORE_dlaset2(PLASMA_enum uplo, int n1, int n2, double alpha,
                     double *tileA, int ldtilea);
    void CORE_dparfb(PLASMA_enum side, PLASMA_enum trans, PLASMA_enum direct, PLASMA_enum storev,
                     int M1, int N1, int M2, int N2, int K, int L,
                           double *A1, int LDA1,
                           double *A2, int LDA2,
                     const double *V, int LDV,
                     const double *T, int LDT,
                           double *WORK, int LDWORK);
}

#endif
