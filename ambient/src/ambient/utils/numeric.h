/*
*copyright Timothee Ewart - University of Geneva, 
*/

#ifndef AMBIENT_UTILS_NUMERIC
#define AMBIENT_UTILS_NUMERIC

#include <alps/hdf5.hpp>

/*
*  idea : partial specialization on the type : double or "complex", complex must be std::complex<double> compatible.
*  The partial specialization adds a safety on the type, function can not be instantiated with float or std::complex<float>
*/

namespace ambient { namespace numeric { namespace kernels {

    template<class T, class D, typename type = void>
    struct helper_cast;

    template<class T, class D>
    struct helper_cast<T, D, typename boost::enable_if< boost::mpl::and_<boost::is_complex<T>, boost::is_floating_point<D> > >::type>{
        static T cast(D a){ return static_cast<T>(a); }
    };

    template<class T, class D>
    struct helper_cast<T, D, typename boost::enable_if< boost::mpl::and_<boost::is_floating_point<T>, boost::is_complex<D> > >::type>{
        static T cast(D a){ return a.real(); }
    };

    template<class T, class D>
    struct helper_cast<T, D, typename boost::enable_if< boost::mpl::and_<boost::is_floating_point<T>, boost::is_floating_point<D> > >::type>{
        static T cast(T a){ return a; }
    };

    template<class T, class D>
    struct helper_cast<T, D, typename boost::enable_if< boost::mpl::and_<boost::is_complex<T>, boost::is_complex<D> > >::type>{
        static T cast(T a){ return a; }
    };

    template <class T>
    inline int OptimalSize(T a){ return (int)a; }

    template <class T>
    inline int OptimalSize(std::complex<T> a){ return (int)a.real(); }

    template<class T, typename type = void>
    struct helper_complex;

    template<class T>
    struct helper_complex<T, typename boost::enable_if<boost::is_floating_point<T> >::type>{
        static inline T conj(T x){ return x; }
        static inline T real(T x){ return x; }
    };
   
    template<class T>
    struct helper_complex<T, typename boost::enable_if<boost::is_complex<T> >::type>{
        static inline T conj(T x){ return std::conj(x); }
        static inline T real(T x){ return std::real(x); }
    };

    template<class t, typename type = void>
    struct helper_blas;

    template<class T>
    struct helper_blas<T, typename boost::enable_if<boost::is_same<T, double> >::type>{ //safety sgemm != dgemm
       static void gemm(const char* transa ,const char* transb, const int* m, const int* n, const int* k, const double* alpha, const T* ad, const int* lda, const T* bd, const int* ldb, const double* beta, T* cd, const int* ldc){
          dgemm_(transa, transb, m, n, k, alpha, ad, lda, bd, ldb, beta, cd, ldc);
       } 
    
       static void gemv(const char* transa, const int* m, const int* n, const double* alfa, const T* ad, const int* lda, const T* bd, const int* ldb, const double* beta, T* cd, const int* ldc){
          dgemv_(transa, m, n, alfa, ad, lda, bd, ldb, beta, cd, ldc);
       }

       static void axpy(const int* n, const double* alpha, const T* ad, const int* inca, T* bd, const int* incb){
          daxpy_(n, alpha, ad, inca, bd, incb);
       }

    };

    template<class T>
    struct helper_blas<T, typename boost::enable_if<boost::mpl::and_<boost::is_complex<T>,boost::is_same<typename T::value_type, double> > >::type>{ //safety cgemm !=  zgemm, your class T must be std::complex<double> compatible
       static void gemm(const char* transa ,const char* transb, const int* m, const int* n, const int* k, const double* alpha, const T* ad, const int* lda, const T* bd, const int* ldb, const double* beta, T* cd, const int* ldc){
          zgemm_(transa, transb, m, n, k, alpha, ad, lda, bd, ldb, beta, cd, ldc);
       } 

       static void gemv(const char* transa, const int* m, const int* n, const double* alfa, const T* ad, const int* lda, const T* bd, const int* ldb, const double* beta, T* cd, const int* ldc){
          zgemv_(transa, m, n, alfa, ad, lda, bd, ldb, beta, cd, ldc);
       }

       static void axpy(const int* n, const T* alfa, const T* ad, const int* inca, T* bd, const int* incb){ // use in exp algo only
          zaxpy_(n, alfa, ad, inca, bd, incb); 
       }

       static void axpy(const int* n, const typename T::value_type* alfa, const T* ad, const int* inca, T* bd, const int* incb){ // use after SVD complex
          T beta = static_cast<T>(*alfa); //double to complex with beta.imag() == 0.0 
          zaxpy_(n, &beta, ad, inca, bd, incb);
       }
    };

    template<class t, typename type = void>
    struct helper_lapack;

    template<class T>
    struct helper_lapack<T, typename boost::enable_if<boost::is_same<T, double> >::type>{
       static void gesvd(const char* jobu, const char* jobvt, const int* m, const int* n, T* ad, const int* lda, T* sd, T* ud, const int* ldu, T* vtd, const int* ldvt, T* wkopt, int* lwork, int* info){
           T* work;
           dgesvd_(jobu, jobvt, m, n, ad, lda, sd, ud, ldu, vtd, ldvt, wkopt, lwork, info);
           *lwork = OptimalSize(*wkopt);
           work = (T*)malloc( (*lwork)*sizeof(T) );
           dgesvd_(jobu, jobvt, m, n, ad, lda, sd, ud, ldu, vtd, ldvt, work, lwork, info);
           assert( *info == 0 ); // otherwise the algorithm computing atomic SVD failed to converge
           free(work);
       } 

       static void syev(const char* jobz, const char* uplo, const int* n, T* a, const int* lda, T* w, T* wkopt, int* lwork, int* info ){
           T* work;
           dsyev_(jobz, uplo, n, a, lda, w, wkopt, lwork, info);
           assert( *info == 0 );
           *lwork = OptimalSize(*wkopt);
           work = (T*)malloc( (*lwork)*sizeof(T) );
           dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
           assert( *info == 0 );
           free(work);
       }

       static void getrf(const int* m, const int*n, T* a, const int* lda, int* ipiv, int* info){
           dgetrf_(m, n, a, lda, ipiv, info);
           assert( *info == 0 );
       }

       static void getri(const int*n, T* a, const int* lda, int* ipiv, int* info){
           T* work;
           T wkopt;
           int lwork = -1;
           dgetri_( n, a, lda, ipiv, &wkopt, &lwork, info);
           assert( *info == 0 );
           lwork = OptimalSize(wkopt);
           work = (T*)malloc(lwork*sizeof(T));
           dgetri_( n, a, lda, ipiv, work, &lwork, info);
           assert( *info == 0 );
           free(work);
       }

       static void larfg(const int *n, T* alpha, T* x, int *incx, T* tau){
           dlarfg_(n, alpha, x, incx, tau); 
           assert( *info == 0 );
       }

       static void gebd2(const int* m, const int* n, T* a, const int* lda, T* d, T* e, T* tauq, T* taup, T* work, int* info){
           dgebd2_(m ,n, a, lda, d, e, tauq, taup, work, info);
           assert( *info == 0 );
       }

       static void gebrd(const int* m, const int* n, T* a, const int* lda, T* d, T* e, T* tauq, T* taup, T* work, const int* lwork, int* info){
           dgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);  
           assert( *info == 0 );
       }

       static void orgbr(const char* vect, const int* m, const int* n, const int* k, T* a, const int* lda, const T* tau, T* work, const int* lwork, int* info){
           dorgbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
           assert( *info == 0 );
       }

       static void gbbrd (const char* vect, const int* m, const int* n, const int* ncc, const int* kl, const int* ku, T* ab, const int* ldab, T* d, T* e, T* q, const int* ldq, T* pt, const int* ldpt, T* c, const int* ldc, T* work, int* info ){
           dggbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
           assert( *info == 0 );
       }
 
       static void bdsqr(const char* uplo, const int* n, const int* ncvt, const int* nru, const int* ncc, T* d, T* e, T* vt, const int* ldvt, T* u, const int* ldu, T* c, const int* ldc, T* work, int* info ){
           dbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info);
           assert( *info == 0 );
       }
    };

    template<class T>
    struct helper_lapack<T, typename boost::enable_if<boost::mpl::and_<boost::is_complex<T>,boost::is_same<typename T::value_type, double> > >::type>{ //safety cgemm !=  zgemm, your class T must be std::complex<double> compatible
       static void gesvd(const char* jobu, const char* jobvt, const int* m, const int* n, T* ad, const int* lda, typename T::value_type* sd, T* ud, const int* ldu, T* vtd, const int* ldvt, T* wkopt, int* lwork, int* info){
           typename T::value_type* rwork = new typename T::value_type[std::max(1,5*std::min((*n),(*m)))]; //from lapack doc
           T* work;
           zgesvd_(jobu, jobvt, m, n, ad, lda, sd, ud, ldu, vtd, ldvt, wkopt, lwork, rwork, info); // query and allocate the optimal workspace
           assert( *info == 0 );
           *lwork = OptimalSize(*wkopt);
           work = (T*)malloc( (*lwork)*sizeof(T) );
           zgesvd_(jobu, jobvt, m, n, ad, lda, sd, ud, ldu, vtd, ldvt, work, lwork, rwork, info); //run
           assert( *info == 0 ); // otherwise the algorithm computing atomic SVD failed to converge
           free(work);
           delete [] rwork;
       } 

       static void syev(const char* jobz, const char* uplo, const int* n, T* a, const int* lda, typename T::value_type* w, T* wkopt, int* lwork, int* info ){
           typename T::value_type* rwork = new typename T::value_type[std::max(1,3*(*n)-2)]; // from intel lapack doc
           T* work;
           zheev_(jobz, uplo, n, a, lda, w, wkopt, lwork, rwork, info);
           assert( *info == 0 );
           *lwork = OptimalSize(*wkopt);
           work = (T*)malloc( (*lwork)*sizeof(T) );
           zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
           assert( *info == 0 );
           free(work);
           delete [] rwork;
       }
  
       static void geev(const char* jobvl, const char* jobvr, const int* n, T* a, const int* lda, T* s, T* ldv, const int* ldlv, T* rvd, const int* ldrv, T* wkopt, int* lwork, int* info){
           typename T::value_type* rwork = new typename T::value_type[std::max(1,2*(*n))]; // from lapack doc
           T* work;
           zgeev_(jobvl, jobvr, n, a, lda, s, ldv, ldlv, rvd, ldrv, wkopt, lwork, rwork, info);
           assert( *info == 0 );
           *lwork = OptimalSize(*wkopt);
           work = (T*)malloc( (*lwork)*sizeof(T) );
           zgeev_(jobvl, jobvr, n, a, lda, s, ldv, ldlv, rvd, ldrv, work, lwork, rwork, info);
           assert( *info == 0 );
           free(work);
           delete [] rwork;
       } 

       static void getrf(const int* m, const int*n, T* a, const int* lda, int* ipiv, int* info){
           zgetrf_(m, n, a, lda, ipiv, info);
       }

       static void getri(const int*n, T* a, const int* lda, int* ipiv, int* info){
           T* work;
           T wkopt;
           int lwork = -1;
           zgetri_( n, a, lda, ipiv, &wkopt, &lwork, info);
           assert( *info == 0 );
           lwork = OptimalSize(wkopt);
           work = (T*)malloc(lwork*sizeof(T));
           zgetri_( n, a, lda, ipiv, work, &lwork, info);
           assert( *info == 0 );
           free(work);
       }

       static void larfg(const int *n, T* alpha, T* x, int *incx, T* tau){
           zlarfg_(n, alpha, x, incx, tau); 
       }

       static void gebd2(const int* m, const int* n, T* a, const int* lda, T* d, T* e, T* tauq, T* taup, T* work, int* info){
           zgebd2_(m ,n, a, lda, d, e, tauq, taup, work, info);
           assert( *info == 0 );
       }

       static void gebrd (const int* m, const int* n, T* a, const int* lda, T * d, T* e, T* tauq, T* taup, T* work,
                        //  typename T::value_type * d, typename T::value_type* e, T* tauq, T* taup, T* work, <------ the good one
                          const int* lwork, int* info){
          assert(false); // fix the signature and the mix double/complex for TE
         //  zgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);  
       }

       static void orgbr(const char* vect, const int* m, const int* n, const int* k, T* a, const int* lda, const T* tau, T* work, const int* lwork, int* info){
           zungbr_(vect, m, n, k, a, lda, tau, work, lwork, info); // double signature != complex signature
       }

       static void gbbrd (const char* vect, const int* m, const int* n, const int* ncc, const int* kl, const int* ku, T* ab, const int* ldab, T* d, T* e, T* q, const int* ldq, T* pt, const int* ldpt, T* c, const int* ldc, T* work, int* info ){
           assert(false); // fix I need one more buffer for complex
           //zggbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
       }

       static void bdsqr(const char* uplo, const int* n, const int* ncvt, const int* nru, const int* ncc, T* d, T* e, T* vt, const int* ldvt, T* u, const int* ldu, T* c, const int* ldc, T* work, int* info ){
          assert(false); // fix the signature and the mix double/complex for TE
        //   zdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info);
       }
    };

    template<class t, typename type = void>
    struct helper_plasma;

    template<class T>
    struct helper_plasma<T, typename boost::enable_if<boost::is_same<T, double> >::type>{
        static void geqrt(int m, int n, int ib, T* a, int lda, T* t, int ldt, T* tau, T* work){
            CORE_dgeqrt(m, n, ib, a, lda, t, ldt, tau, work);  
        }

        static void ormqr(int side, int trans, int m, int n, int k, int in, const T *a, int lda, const T *t, int ldt, T *c, int ldc, T *work, int ldwork){
            CORE_dormqr(side, trans, m, n, k, in, a, lda, t, ldt, c, ldc, work, ldwork);
        }
       
        static void tsqrt(int m, int n, int ib, T* a1, int lda1, T* a2, int lda2, T* t, int ldt, T* tau, T* work){
            CORE_dtsqrt(m, n, ib, a1, lda1, a2, lda2, t, ldt, tau, work);
        }

        static void tsmqr(int side, int trans, int m1, int n1, int m2, int n2, int k, int ib, T* a1, int lda1, T* a2, int lda2, T* v, int ldv, T* t, int ldt, T* work, int ldwork){
            CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, a1, lda1, a2, lda2, v, ldv, t, ldt, work, ldwork);
        }

        static void gelqt(int m, int n, int ib, T* a, int lda, T* t, int ldt, T* tau, T* work){
            CORE_dgelqt(m, n, ib, a, lda, t, ldt, tau, work);  
        }

        static void ormlq(int side, int trans, int m, int n, int k, int in, const T *a, int lda, const T *t, int ldt, T *c, int ldc, T *work, int ldwork){
            CORE_dormlq(side, trans, m, n, k, in, a, lda, t, ldt, c, ldc, work, ldwork);
        }

        static void tslqt(int m, int n, int ib, T* a1, int lda1, T* a2, int lda2, T* t, int ldt, T* tau, T* work){
            CORE_dtslqt(m, n, ib, a1, lda1, a2, lda2, t, ldt, tau, work);
        }

        static void tsmlq(int side, int trans, int m1, int n1, int m2, int n2, int k, int ib, T* a1, int lda1, T* a2, int lda2, T* v, int ldv, T* t, int ldt, T* work, int ldwork){
            CORE_dtsmlq(side, trans, m1, n1, m2, n2, k, ib, a1, lda1, a2, lda2, v, ldv, t, ldt, work, ldwork);
        }

        static void laset2(int uplo, int n1, int n2, T alpha, T* A, int lda){
            CORE_dlaset2(uplo, n1, n2, alpha, A, lda);
        }

    };

    template<class T>
    struct helper_plasma<T, typename boost::enable_if<boost::mpl::and_<boost::is_complex<T>,boost::is_same<typename T::value_type, double> > >::type>{ 
        static void geqrt(int m, int n, int ib, T* a, int lda, T* t, int ldt, T* tau, T* work){
            CORE_zgeqrt(m, n, ib, a, lda, t, ldt, tau, work);  
        }

        static void ormqr(int side, int trans, int m, int n, int k, int in, const T *a, int lda, const T *t, int ldt, T *c, int ldc, T *work, int ldwork){
            CORE_zunmqr(side, trans, m, n, k, in, a, lda, t, ldt, c, ldc, work, ldwork);
        }

        static void tsqrt(int m, int n, int ib, T* a1, int lda1, T* a2, int lda2, T* t, int ldt, T* tau, T* work){
            CORE_ztsqrt(m, n, ib, a1, lda1, a2, lda2, t, ldt, tau, work);
        }

        static void tsmqr(int side, int trans, int m1, int n1, int m2, int n2, int k, int ib, T* a1, int lda1, T* a2, int lda2, T* v, int ldv, T* t, int ldt, T* work, int ldwork){
            CORE_ztsmqr(side, trans, m1, n1, m2, n2, k, ib, a1, lda1, a2, lda2, v, ldv, t, ldt, work, ldwork);
        }

        static void gelqt(int m, int n, int ib, T* a, int lda, T* t, int ldt, T* tau, T* work){
            CORE_zgelqt(m, n, ib, a, lda, t, ldt, tau, work);  
        }

        static void ormlq(int side, int trans, int m, int n, int k, int in, const T *a, int lda, const T *t, int ldt, T *c, int ldc, T *work, int ldwork){
            CORE_zunmlq(side, trans, m, n, k, in, a, lda, t, ldt, c, ldc, work, ldwork);
        }

        static void tslqt(int m, int n, int ib, T* a1, int lda1, T* a2, int lda2, T* t, int ldt, T* tau, T* work){
            CORE_ztslqt(m, n, ib, a1, lda1, a2, lda2, t, ldt, tau, work);
        }

        static void tsmlq(int side, int trans, int m1, int n1, int m2, int n2, int k, int ib, T* a1, int lda1, T* a2, int lda2, T* v, int ldv, T* t, int ldt, T* work, int ldwork){
            CORE_ztsmlq(side, trans, m1, n1, m2, n2, k, ib, a1, lda1, a2, lda2, v, ldv, t, ldt, work, ldwork);
        }

        static void laset2(int uplo, int n1, int n2, T alpha, T* A, int lda){
            CORE_zlaset2(uplo, n1, n2, alpha, A, lda);
        }
    };

}}} // end namespace

#endif

