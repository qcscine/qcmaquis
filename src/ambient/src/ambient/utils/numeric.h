#ifndef NUMERIC_H
#define NUMERIC_H

#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/control/expr_if.hpp>

// C - dgemm, zgemm = prefix + gemm so tuples !
#define ARITHMETIC_TYPE(I) ARITHMETIC_TYPE ## I
#define TUPLE_TYPE0 (double,d)
#define TUPLE_TYPE1 (void,z)
// C - the wrapper checks the type
#define TUPLE_WRAPPER0 (double,d)
#define TUPLE_WRAPPER1 (std::complex<double>,z)

#define TUPLE_TYPE_CNT 2 

// C - Lapack complex solver e.g. zgesvd needs one more arguments
#define COMMA_(I) BOOST_PP_IF(I, BOOST_PP_COMMA, BOOST_PP_EMPTY)()
#define ARG(I) BOOST_PP_EXPR_IF(I, double* rwork)

// C - Two functions for the validation tests
    inline double norm(double a){return sqrt(a*a);}
    inline double norm(std::complex<double> a){return sqrt(std::norm(a));}

// C - Two functions for the SVD, to calculate the optimal size work because a cast complex to int does not exist !
    inline int OptimalSize(double a){return (int)a;}
    inline int OptimalSize(std::complex<double> a){return (int)a.real();}

extern "C" {
// C - it should be useless, to do find the bug !
    double sqrt(double);
    double fabs(double);

// C - Fortran declaration : native compatibility all math libs (MKL, Goto, ACML, ESSLSMP, ...)
#define BOOST_PP_DEF(z, I, _) \
    void BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,1,TUPLE_TYPE ## I),gemm_)(const char *transa, const char *transb, const int  *m, const int *n, const int *k, \
                      const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *alpha, const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *a, const int *lda, const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *b, const int *ldb, \
                      const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *beta, BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *c, const int *ldc); \
    void BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,1,TUPLE_TYPE ## I),axpy_)(const int *n, const  BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *alpha, const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *x, const int *incx, \
                      BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *y, const int *incy); \
    void BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,1,TUPLE_TYPE ## I),gesvd_)( const char* jobu, const char* jobvt, const int* m, \
                      const int* n, BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I) *a, const int* lda, double* s, \
                      BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I)* u, const int* ldu, BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I)* vt, const int* ldvt, \
                      BOOST_PP_TUPLE_ELEM(2,0,TUPLE_TYPE ## I)* work, const int* lwork, ARG(I) COMMA_(I) int* info );

BOOST_PP_REPEAT(TUPLE_TYPE_CNT, BOOST_PP_DEF, _)
#undef BOOST_PP_DEF

    // C - no complex version for this function
    void dsyev_( const char* jobz, const char* uplo, const int* n, double* a, 
                 const int* lda, double* w, double* work, const int* lwork, 
                 int* info );
}
/*
* wrapper to the lib call 
*/
#define BOOST_PP_DEF(z, I, _) \
    inline void gemm(const char *transa, const char *transb, const int  *m, const int *n, const int *k, \
                     const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *alpha, const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *a, const int *lda, const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *b, const int *ldb, \
                     const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *beta, BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *c, const int *ldc) \
                { \
                    BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,1,TUPLE_WRAPPER ## I),gemm_)(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); \
                } \
    inline void axpy(const int *n, const  BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *alpha, const BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *x, const int *incx, \
                     BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *y, const int *incy) \
                { \
                    BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,1,TUPLE_WRAPPER ## I),axpy_)(n, alpha, x, incx, y, incy); \
                } \
    inline void gesvd(const char* jobu, const char* jobvt, const int* m, \
                      const int* n, BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I) *a, const int* lda, double * s, \
                      BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I)* u, const int* ldu, BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I)* vt, const int* ldvt, \
                      BOOST_PP_TUPLE_ELEM(2,0,TUPLE_WRAPPER ## I)* work, const int* lwork, double* rwork, int* info ) \
                { \
                    BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,1,TUPLE_WRAPPER ## I),gesvd_)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, BOOST_PP_EXPR_IF(I,rwork) COMMA_(I) info); \
                }
    /**/
BOOST_PP_REPEAT(TUPLE_TYPE_CNT, BOOST_PP_DEF, _)
#undef BOOST_PP_DEF

#endif
