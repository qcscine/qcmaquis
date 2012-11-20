#ifndef AMBIENT_NUMERIC_KERNELS_UTILS
#define AMBIENT_NUMERIC_KERNELS_UTILS

#include <limits>
#include "utils/timings.h"

//#define AMBIENT_COMPUTATIONAL_TIMINGS
//#define AMBIENT_CHECK_BOUNDARIES

#ifdef AMBIENT_CHECK_BOUNDARIES
#include <execinfo.h>
#endif
    
#include "ambient/utils/ceil.h"

#define NUM_ROWS(m,offset) std::min(IB, m.num_rows() - offset % m.num_rows())
#define NUM_COLS(m,offset) std::min(IB, m.num_cols() - (size_t)(offset / m.num_rows()))

#ifdef AMBIENT_COMPUTATIONAL_TIMINGS
    #define __A_TIME_C(name) static __a_timer time(name); time.begin();
    #define __A_TIME_C_STOP time.end();
#else
    #define __A_TIME_C(name) 
    #define __A_TIME_C_STOP 
#endif

typedef int PLASMA_enum;
#define PlasmaNoTrans    111
#define PlasmaTrans      112
#define PlasmaUpper      121
#define PlasmaLower      122
#define PlasmaUpperLower 123
#define PlasmaLeft       141
#define PlasmaRight      142
#define PlasmaForward    391
#define PlasmaColumnwise 401

extern "C" {
void      dgemm_(const char*,const char*, const int*, const int*, const int*, const double*, 
                 const double*, const int*, const double*, const int*, 
                 const double*, double*, const int*);
void     dbdsqr_(const char*, const int*, const int*, const int*, const int*, double*, double*, 
                 double*, const int*, double*, const int*, double*, const int*, double*, int*);

int LAPACKE_dbdsdc(int matrix_order, char uplo, char compq, int n, double* d, double* e, double* u, int ldu, double* vt, int ldvt, double* q, int* iq);

void     dgebrd_(const int*, const int*, double*, const int*, double*, double*, double*, double*, double*, const int*, int*);
double     ddot_(const int*, const double*, const int*, 
                 const double*, const int*);
void     dgebd2_(int*, int*, double*, int*, double*, double*, double*, double*, double*, int*);
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

namespace ambient { namespace numeric { namespace kernels {

    #include "ambient/utils/numeric.h"

    template<typename T> inline dim2 __a_get_dim(T& ref){ 
        return ref.impl->spec.dim;
    }

    template<typename T> inline size_t __a_sizeof(T& ref){ 
        return ref.impl->spec.size;  
    }

    inline void* __a_solidify(char* r, size_t sz){
        void* memory = malloc(sz);
        memcpy(memory, r, sz);
        return memory;
    }

    inline void __a_disperse(void* data, char* r, size_t sz){
        memcpy(r, data, sz);
        free(data);
    }

    template <typename T> inline T __a_dot(T* a, T* b, int size){
        T summ(0);
        for(size_t k=0; k < size; k++)
           summ += a[k]*b[k];
        return summ;
    }

    inline double __a_dot(double* a, double* b, int size){
        static const int ONE = 1;
        return ddot_(&size, a, &ONE, b, &ONE);
    }

    template <typename T>
    inline void __a_memcpy(T* dd, T *sd, size_t w, T alfa){
        memcpy(dd, sd, w);
    }

    template <typename T>
    inline void __a_memscal(T* dd, T *sd, size_t w, T alfa){
        int z = w/sizeof(T);
        do{ *dd++ = alfa*(*sd++); }while(--z > 0); // be carefull that dd != sd
    }

    template <typename T>
    inline void __a_memscala(T* dd, T *sd, size_t w, T alfa){
        int z = w/sizeof(T);
        do{ *dd++ += alfa*(*sd++); }while(--z > 0); // be carefull that dd != sd
    }

    template<typename T, void(*PTF)(T* dd, T* sd, size_t w, T alfa)>
    inline void __a_memptf_r(T* dst, int ldb, dim2 dst_p, 
                                    T* src, int lda, dim2 src_p, 
                                    dim2 size, T alfa = 0.0)
    {
#ifdef AMBIENT_CHECK_BOUNDARIES
        if(__a_get_dim(dst).x - dst_p.x < size.x || __a_get_dim(dst).y - dst_p.y < size.y ||
           __a_get_dim(src).x - src_p.x < size.x || __a_get_dim(src).y - src_p.y < size.y){
            ambient::cout << "Error: invalid memory movement: " << std::endl;
            ambient::cout << "Matrix dst " << __a_get_dim(dst).x << "x" << __a_get_dim(dst).y << "\n";
            ambient::cout << "Dest p " << dst_p.x << "x" << dst_p.y << "\n";
            ambient::cout << "Matrix src " << __a_get_dim(src).x << "x" << __a_get_dim(src).y << "\n";
            ambient::cout << "Src p " << src_p.x << "x" << src_p.y << "\n";
            ambient::cout << "Block size " << size.x << "x" << size.y << "\n";

            void *array[10];
            size_t size = backtrace(array, 10);
            backtrace_symbols_fd(array, size, 2);
        }
#endif
        int n = size.x;
        int m = size.y*sizeof(T);

        T* sd = src + src_p.y + src_p.x*lda;
        T* dd = dst + dst_p.y + dst_p.x*ldb;

        do{ PTF(dd, sd, m, alfa); sd += lda; dd += ldb; }while(--n > 0);
    }

    template <typename T>
    inline void __a_refresh(T* dst, const T* src, size_t sz){
        if(dst != src) memcpy(dst, src, sz);
    }

} } }

#endif
