/*
*copyright Timothee Ewart - University of Geneva, 
*/

#ifndef AMBIENT_UTILS_NUMERIC
#define AMBIENT_UTILS_NUMERIC

/*
*  idea : partial specialization on the type : double or "complex", complex must be std::complex<double> compatible.
*  The partial specialization adds a safety on the type, function can not be instantiated with float or std::complex<float>
*/

namespace ambient { namespace numeric { namespace kernels {

    inline int OptimalSize(double a){return (int)a;}
    inline int OptimalSize(std::complex<double> a){return (int)a.real();}

    template<class t, typename type = void>
    struct helper_blas;

    template<class T>
    struct helper_blas<T, typename boost::enable_if<boost::is_same<T, double> >::type>{ //safety sgemm != dgemm
       static void gemm(const char* transa ,const char* transb, const int* m, const int* n, const int* k,
                        const double* alpha, const T* ad, const int* lda, const T* bd, const int* ldb, const double* beta, T* cd, const int* ldc){
          dgemm_(transa, transb, m, n, k, alpha, ad, lda, bd, ldb, beta, cd, ldc);
       } 
    
       static void gemv(const char* transa, const int* m, const int* n,
                        const double* alpha, const T* ad, const int* lda, const T* bd, const int* ldb, const double* beta, const T* cd, const int* ldc){
          dgemv_(transa, m, n, alfa, ad, lda, bd, ldb, beta, cd, ldc);
       }

       static void axpy(const int* n, const double* alpha, const T* ad, const int* inca, T* bd, const int* incb){
          daxpy_(n, alpha, ad, inca, bd, incb);
       }

    };

    template<class T>
    struct helper_blas<T, typename boost::enable_if<boost::mpl::and_<boost::is_complex<T>,boost::is_same<typename T::value_type, double> > >::type>{ //safety cgemm !=  zgemm, your class T must be std::complex<double> compatible
       static void gemm(const char* transa ,const char* transb, const int* m, const int* n, const int* k,
                        const double* alpha, const T* ad, const int* lda, const T* bd, const int* ldb, const double* beta, T* cd, const int* ldc){
          zgemm_(transa, transb, m, n, k, alpha, ad, lda, bd, ldb, beta, cd, ldc);
       } 

       static void gemv(const char* transa, const int* m, const int* n,
                        const double* alpha, const T* ad, const int* lda, const T* bd, const int* ldb, const double* beta, const T* cd, const int* ldc){
          zgemv_(transa, m, n, alfa, ad, lda, bd, ldb, beta, cd, ldc);
       }

       static void axpy(const int* n, const T* alpha, const T* ad, const int* inca, T* bd, const int* incb){
          zaxpy_(n, alpha, ad, inca, bd, incb);
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
           assert( info == 0 ); // otherwise the algorithm computing atomic SVD failed to converge
           free(work);
       } 

       static void syev(const char* jobz, const char* uplo, const int* n, T* a, const int* lda, T* w, T* wkopt, int* lwork, int* info ){
           T* work;
           dsyev_(jobz, uplo, n, a, lda, w, wkopt, lwork, info);
           *lwork = OptimalSize(*wkopt);
           work = (T*)malloc( (*lwork)*sizeof(T) );
           dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
           assert( info == 0 );
           free(work);
       }

       static void larfg(const int *n, T* alpha, T* x, int *incx, T* tau){
           dlarfg_(n, alpha, x, incx, tau); 
       }

       static void gebd2(const int* m, const int* n, T* a,
                         const int* lda, T* d, T* e, T* tauq,
                         T* taup, T* work, int* info){
           dgebd2_(m ,n, a, lda, d, e, tauq, taup, work, info);
       }

       static void gebrd(const int* m, const int* n, T* a, const int* lda,
                          T* d, T* e, T* tauq, T* taup, T* work,
                          const int* lwork, int* info){
           dgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);  
       }

       static void orgbr(const char* vect, const int* m, const int* n, 
                          const int* k, T* a, const int* lda,
                          const T* tau, T* work, const int* lwork,
                          int* info){
           dorgbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
       }

       static void gbbrd (const char* vect, const int* m, const int* n, 
                          const int* ncc, const int* kl, const int* ku,
                          T* ab, const int* ldab, T* d, T* e, T* q,
                          const int* ldq, T* pt, const int* ldpt, T* c,
                          const int* ldc, T* work, int* info ){
           dggbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
       }
 
       static void bdsqr(const char* uplo, const int* n, const int* ncvt,
                         const int* nru, const int* ncc, T* d, T* e,
                         T* vt, const int* ldvt, T* u, const int* ldu,
                         T* c, const int* ldc, T* work, int* info ){
           dbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info);
       }
    };

    template<class T>
    struct helper_lapack<T, typename boost::enable_if<boost::mpl::and_<boost::is_complex<T>,boost::is_same<typename T::value_type, double> > >::type>{ //safety cgemm !=  zgemm, your class T must be std::complex<double> compatible
       static void gesvd(const char* jobu, const char* jobvt, const int* m, const int* n, T* ad, const int* lda, typename T::value_type* sd, T* ud, const int* ldu, T* vtd, const int* ldvt, T* wkopt, int* lwork, int* info){
           typename T::value_type* rwork = new typename T::value_type[5*(*m)]; //from lapack doc
           T* work;
           zgesvd_(jobu, jobvt, m, n, ad, lda, sd, ud, ldu, vtd, ldvt, wkopt, lwork, rwork, info); // query and allocate the optimal workspace
           *lwork = OptimalSize(*wkopt);
           work = (T*)malloc( (*lwork)*sizeof(T) );
           zgesvd_(jobu, jobvt, m, n, ad, lda, sd, ud, ldu, vtd, ldvt, work, lwork, rwork, info); //run
           assert( info == 0 ); // otherwise the algorithm computing atomic SVD failed to converge
           free(work);
           delete [] rwork;
       } 

       static void syev(const char* jobz, const char* uplo, const int* n, T* a, const int* lda, typename T::value_type* w, T* wkopt, int* lwork, int* info ){
           typename T::value_type* rwork = new typename T::value_type[std::max(1,2*(*n))]; // from lapack doc
           T* work;
           zheev_(jobz, uplo, n, a, lda, w, wkopt, lwork, rwork, info);
           *lwork = OptimalSize(*wkopt);
           work = (T*)malloc( (*lwork)*sizeof(T) );
           zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
           assert( info == 0 );
           free(work);
           delete [] rwork;
       };

       static void larfg(const int *n, T* alpha, T* x, int *incx, T* tau){
           zlarfg_(n, alpha, x, incx, tau); 
       }

       static void gebd2(const int* m, const int* n, T* a,
                         const int* lda, T* d, T* e, T* tauq,
                         T* taup, T* work, int* info){
           zgebd2_(m ,n, a, lda, d, e, tauq, taup, work, info);
       }

       static void gebrd (const int* m, const int* n, T* a, const int* lda,
                          T * d, T* e, T* tauq, T* taup, T* work,
                        //  typename T::value_type * d, typename T::value_type* e, T* tauq, T* taup, T* work, <------ the good one
                          const int* lwork, int* info){
          assert(false); // fix the signature and the mix double/complex for TE
         //  zgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);  
       };

       static void orgbr(const char* vect, const int* m, const int* n, 
                          const int* k, T* a, const int* lda,
                          const T* tau, T* work, const int* lwork,
                          int* info){
           zungbr_(vect, m, n, k, a, lda, tau, work, lwork, info); // double signature != complex signature
       };

       static void gbbrd (const char* vect, const int* m, const int* n, 
                          const int* ncc, const int* kl, const int* ku,
                          T* ab, const int* ldab, T* d, T* e, T* q,
                          const int* ldq, T* pt, const int* ldpt, T* c,
                          const int* ldc, T* work, int* info ){
           assert(false); // fix I need one more buffer for complex
           //zggbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
       }

       static void bdsqr(const char* uplo, const int* n, const int* ncvt,
                         const int* nru, const int* ncc, T* d, T* e,
                         T* vt, const int* ldvt, T* u, const int* ldu,
                         T* c, const int* ldc, T* work, int* info ){
          assert(false); // fix the signature and the mix double/complex for TE
        //   zdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info);
       }
    };
}}} // end namespace

#endif

