#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "ambient/numeric/matrix.hpp"
#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "utilities.h"

extern "C" {
    void dgemv(const char* trans, int* m, int* n, double* alfa, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
    void dlarfg(int* n, double* alfa, double* x, int* incx, double* tau);
    void dscal(int* n, double* da, double* dx, int* incx);
}

void mkl_labrd( int m, int n, int nb, double* a, int lda, 
            double* d, double* e, double* tauq, double* taup, 
            double* x, int ldx, double* y, int ldy )
{
    double mone = -1.;
    double one = 1.;
    double zero = 0.;
    int lone = 1;
    if(m >= n){ // Reduce to upper bidiagonal form
        for(int i = 0; i < nb; ++i){

            // Update a(i:m,i)
            int ri  = m-i;   //std::min(m-i, i*nb);
            int rj  = n-i-1; //std::min(n-i-1, (i+1)*nb);
            int rij = m-i-1; //std::min(m-i-1, (i+1)*nb);
            dgemv_("N", &ri, &i, &mone, &a[ i ], &lda, &y[ i ], &ldy, &one, &a[i + i*lda], &lone);
            dgemv_("N", &ri, &i, &mone, &x[ i ], &ldx, &a[ i*lda ], &lone, &one, &a[i + i*lda], &lone);

            // Generate reflection Q(i) to annihilate a(i+1:m,i)
            dlarfg_( &ri, &a[i+i*lda], &a[std::min(i+1, m-1)+i*lda], &lone, &tauq[i] );
            d[i] = a[i+i*lda];

            if(i < n-1){
                a[i+i*lda] = one;

                // Compute y(i+1:n,i)
                dgemv_("T", &ri, &rj, &one, &a[i + (i+1)*lda], &lda, &a[i + i*lda], &lone, &zero, &y[i+1 + i*ldy], &lone );
                dgemv_("T", &ri, &i, &one, &a[ i ], &lda, &a[i +i*lda], &lone, &zero, &y[i*ldy], &lone);
                dgemv_("N", &rj, &i, &mone, &y[ i+1 ], &ldy, &y[i*ldy], &lone, &one, &y[i+1 + i*ldy], &lone);
                dgemv_("T", &ri, &i, &one, &x[i], &ldx, &a[i + i*lda], &lone, &zero, &y[ i*ldy ], &lone);
                dgemv_("T", &i, &rj, &mone, &a[ (i+1)*lda ], &lda, &y[i*ldy], &lone, &one, &y[ i+1 + i*ldy], &lone);


                //     n-i-1 ::    y(i+1,i) * tau(i)
                dscal_(&rj, &tauq[i], &y[i+1 + i*ldy], &lone);
                y[i+i*ldy] = zero;

                // Update a(i,i+1:n)
                int r3 = i+1;
                dgemv_("N", &rj, &r3, &mone, &y[ i+1 ], &ldy, &a[i], &lda, &one, &a[i + (i+1)*lda], &lda);
                dgemv_("T", &i, &rj, &mone, &a[(i+1)*lda], &lda, &x[i], &ldx, &one, &a[ i + (i+1)*lda], &lda);

                // Generate reflection P(i) to annihilate a(i,i+2:n)
                dlarfg_(&rj, &a[i + (i+1)*lda], &a[i + std::min( i+2, n-1 )*lda], &lda, &taup[i] );
                e[i] = a[i + (i+1)*lda];
                a[i + (i+1)*lda] = one;

                // Compute X(i+1:m,i)
                dgemv_("N", &rij, &rj, &one, &a[i+1+(i+1)*lda], &lda, &a[i +(i+1)*lda], &lda, &zero, &x[i+1+i*ldx], &lone);
                dgemv_("T", &rj, &r3, &one, &y[i+1], &ldy, &a[i+(i+1)*lda], &lda, &zero, &x[i*ldx], &lone);
                dgemv_("N", &rij, &r3, &mone, &a[i+1], &lda, &x[i*ldx], &lone, &one, &x[i+1+i*ldx], &lone);
                dgemv_("N", &i, &rj, &one, &a[(i+1)*lda], &lda, &a[ i +(i+1)*lda], &lda, &zero, &x[i*ldx], &lone);
                dgemv_("N", &rij, &i, &mone, &x[i+1], &ldx, &x[i*ldx], &lone, &one, &x[i+1+i*ldx], &lone);

                dscal_(&rij, &taup[i], &x[i+1+ i*ldx], &lone );
            }
        }
    }else{}
}

void mkl_gebrd(int m, int n, int nb, double* a, int lda, double* d, double* e, double* tauq, double* taup, double* work, int lwork){
    int minmn = std::min(m,n);
    int ldwrkx = m;
    int ldwrky = n;
    int nx = minmn;
    double one = 1.;
    double mone = -1.;
    int info;
    int i;

    if(nb < minmn){ nx %= nb; if(nx == 0) nx = nb; }
    for(i = 0; i < minmn - nx; i += nb){
        int tm = m-i-nb;
        int tn = n-i-nb;
        mkl_labrd(m-i, n-i, nb, &a[i+i*lda], lda, &d[i], &e[i], &tauq[i], &taup[i], work, ldwrkx, &work[ldwrkx*nb], ldwrky);
        dgemm_( "N", "T", &tm, &tn, &nb, &mone, &a[i+nb+i*lda], &lda, &work[ldwrkx*nb+nb], &ldwrky, &one, &a[i+nb+(i+nb)*lda], &lda);
        dgemm_( "N", "N", &tm, &tn, &nb, &mone, &work[nb], &ldwrkx, &a[i+(i+nb)*lda], &lda, &one, &a[i+nb+(i+nb)*lda], &lda);
    }
    int um = m-i;
    int un = n-i;
    dgebd2_(&um, &un, &a[i+i*lda], &lda, &d[i], &e[i], &tauq[i], &taup[i], work, &info);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( BIDIAGONALIZATION, T, test_types)
{
    pMatrix A(T::valuex,T::valuey);
    pMatrix U(T::valuex,T::valuey);
    pMatrix V(T::valuex,T::valuey);
    pMatrix S(T::valuex,T::valuey);

    generate(A);
    ambient::numeric::band(A,U,S,V);

#ifdef PGEBRD   
    int nb = AMBIENT_IB; 
    int m = num_rows(S);
    int n = num_cols(S);
    int k = std::min(m,n);

    pDiagMatrix d(k);
    pDiagMatrix e(k); 
    pMatrix u(k,k); 
    pMatrix v(k,k); 

    pgebrd(S, d, e, u, v);

    merge(S);
    printf("\n\n\nAfter banding!\n\n");
    std::cout << S[0] << "\n";

    merge(d);
    printf("\n\nD:\n");
    std::cout << d[0] << "\n";

    merge(e);
    printf("\n\nE:\n");
    std::cout << e[0] << "\n";
#elif defined GEBRD
    merge(S);

    int nb = 4;
    int m = S.num_rows();
    int n = S.num_cols();
    int k = std::min(m,n);
    int lda = S.num_rows();
    int lwork = (m+n)*nb;
    double* a = &S[0](0,0); 
    double* d    = (double*)malloc(sizeof(double)*k);     memset(d,    0, sizeof(double)*k);
    double* e    = (double*)malloc(sizeof(double)*k);     memset(e,    0, sizeof(double)*k);
    double* taup = (double*)malloc(sizeof(double)*k);     memset(taup, 0, sizeof(double)*k);
    double* tauq = (double*)malloc(sizeof(double)*k);     memset(tauq, 0, sizeof(double)*k);
    double* work = (double*)malloc(sizeof(double)*lwork); memset(work, 0, sizeof(double)*lwork);

    mkl_gebrd( m, n, nb, a, lda, d, e, tauq, taup, work, lwork);

    printf("\n\n\nAfter banding!\n\n");
    std::cout << S[0] << "\n";

    printf("\n\nTAUP:\n");
    for(int i = 0; i < k; i++)
        printf("%.2f \n\n", taup[i]);

    printf("\n\nTAUQ:\n");
    for(int i = 0; i < k; i++)
        printf("%.2f \n\n", tauq[i]);

    printf("\n\nD:\n");
    for(int j = 0; j < k; j++)
        printf("%.2f \n\n", d[j]);

    printf("\n\nE:\n");
    for(int j = 0; j < k; j++)
        printf("%.2f \n\n", e[j]);
#endif

    pMatrix G1(num_rows(U),num_cols(S));
    pMatrix G2(num_rows(U),num_cols(V));

    ambient::numeric::gemm(U, S, G1);
    ambient::numeric::gemm(G1, V, G2);

    std::cout << "Done: " << num_rows(A) << " " << num_cols(A) << "\n";
    BOOST_CHECK(G2 == A);
}
