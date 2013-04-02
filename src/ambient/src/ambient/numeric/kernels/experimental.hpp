#ifndef AMBIENT_NUMERIC_MATRIX_KERNELS_EXPERIMENTAL
#define AMBIENT_NUMERIC_MATRIX_KERNELS_EXPERIMENTAL

#include "ambient/numeric/kernels/math.hpp"
#include "ambient/numeric/kernels/utils.hpp"

namespace ambient { namespace numeric { namespace kernels {


    template<typename T, PLASMA_enum UL, size_t OFF>
    struct laset2 : public kernel< laset2<T,UL,OFF> > 
    { static void c(matrix<T>& a, const double& alfa); };

    template<int alfa, typename T>
    struct add_vectors : public kernel< add_vectors<alfa, T> > 
    { static void c(matrix<T>& a, const size_t& aoffset, const matrix<T>& b, const size_t& boffset, const size_t& size); };

    template<typename T>
    struct labrd_update_col : public kernel< labrd_update_col<T> > 
    { static void c(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, matrix<T>& tq, matrix<T>& d, const int& i); };

    template<typename T>
    struct labrd_reduce_col : public kernel< labrd_reduce_col<T> > 
    { static void c(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, const int& i); };

    template<typename T>
    struct labrd_update_row : public kernel< labrd_update_row<T> > 
    { static void c(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, matrix<T>& tp, matrix<T>& e, const int& i); };

    template<typename T>
    struct labrd_reduce_row : public kernel< labrd_reduce_row<T> > 
    { static void c(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, const int& i); };

    template<typename T, PLASMA_enum TR>
    struct larfg : public kernel< larfg<T,TR> > 
    { static void c(matrix<T>& a, matrix<T>& t, matrix<T>& d, const size_t& k); };

    template<typename T>
    struct gebd2 : public kernel< gebd2<T> > 
    { static void c(matrix<T>& a, matrix<T>& d, matrix<T>& e, matrix<T>& tq, matrix<T>& tp); };

    template<typename T>
    struct gebrd : public kernel< gebrd<T> > 
    { static void c(matrix<T>& a, weak_view<T>& d, weak_view<T>& e, weak_view<T>& q, weak_view<T>& p); };

    template<typename T>
    struct gbbrd : public kernel< gbbrd<T> > 
    { static void c(matrix<T>& a, weak_view<T>& d, weak_view<T>& e, weak_view<T>& q, weak_view<T>& p); };

    template<typename T>
    struct bdsqr : public kernel< bdsqr<T> > 
    { static void c(matrix<T>& d, matrix<T>& e, matrix<T>& u, matrix<T>& v); };

    template<typename T, PLASMA_enum UL>
    struct copy_band : public kernel< copy_band<T,UL> > 
    { static void c(const matrix<T>& src, matrix<T>& dst, const size_t& dj); };

    template<int ADD, class VA, class VB, class VC, class VF>
    struct gemv_scale : public kernel< gemv_scale<ADD, VA, VB, VC, VF> > {
        typedef typename VC::value_type T;
        static void c(const matrix<T>& a, const size_t& aoffset, 
                      const matrix<T>& b, const size_t& boffset,
                            matrix<T>& c, const size_t& coffset,
                      const matrix<T>& f, const size_t& foffset,
                      const size_t& rows, const size_t& cols);
    };

    template<int alfa, int beta, class ViewA, class ViewB, class ViewC>
    struct gemv : public kernel< gemv<alfa, beta, ViewA, ViewB, ViewC> > {
        typedef typename ViewC::value_type T;
        static void c(const matrix<T>& a, const size_t& aoffset, 
                      const matrix<T>& b, const size_t& boffset,
                            matrix<T>& c, const size_t& coffset,
                      const size_t& rows, const size_t& cols);
    };




/////////////////////////////////////////
// Experimental kernels implementation //
/////////////////////////////////////////



    template<typename T, PLASMA_enum UL, size_t OFF>
    void laset2<T,UL,OFF>::c(matrix<T>& a, const double& alfa){
        double* ad = (double*)s_updated(a);
        CORE_dlaset2(UL, a.num_rows()-OFF, a.num_cols()-OFF, alfa, ad + OFF*a.num_rows(), a.num_rows());
    }

    template<int alfa, typename T>
    void add_vectors<alfa, T>::c(matrix<T>& a, const size_t& aoffset, const matrix<T>& b, const size_t& boffset, const size_t& size){
        T* ad = &((T*)current(a))[aoffset];
        T* bd = &((T*)current(b))[boffset];
        T* ar = &((T*)updated(a))[aoffset];

        for(int k = 0; k < size; k++) 
            ar[k] = alfa*ad[k] + bd[k];
    }
        
    template<typename T>
    void labrd_update_col<T>::c(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, matrix<T>& tq, matrix<T>& d, const int& i){
        static const double mone = -1.;
        static const double one = 1.;
        static const double zero = 0.;
        static const int lone = 1;

        int m  = num_rows(say);
        int n  = num_cols(sax);
        int ri = m-i;
        int rj = n-i-1;

        T* sayd = s_updated(say); int ldsay = say.num_rows();
        T* saxd = current(sax);   int ldsax = sax.num_rows();
        T* syd  = s_updated(sy);  int ldsy = sy.num_rows();
        T* sxd  = current(sx);    int ldsx = sx.num_rows();
        T* tqd  = s_updated(tq);
        T* dd   = s_updated(d);
        
        if(i == 0){
            dlarfg_(&ri, sayd, &sayd[1], &lone, tqd);
            *dd = *sayd;
            *sayd = 1.0;
            return;
        }

        ambient::memptf<T, ambient::memcpy>(sayd, ldsay, dim2(i, i-1), 
                                            saxd, ldsax, dim2(i, i-1), 
                                            dim2( num_cols(say)-i, 1));

        dgemv_("N", &ri, &i, &mone, &sayd[ i ], &ldsay, &syd[ i ], &ldsy, &one, &sayd[i + i*ldsay], &lone);
        dgemv_("N", &ri, &i, &mone, &sxd[ i ], &ldsx, &sayd[ i*ldsay ], &lone, &one, &sayd[i + i*ldsay], &lone);
        
        dlarfg_( &ri, &sayd[i+i*ldsay], &sayd[std::min(i+1, m-1)+i*ldsay], &lone, &tqd[i] );
        dd[i] = sayd[i+i*ldsay];
        sayd[i+i*ldsay] = 1.000;
    }

    template<typename T>
    void labrd_reduce_col<T>::c(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, const int& i){
        static const double mone = -1.;
        static const double one = 1.;
        static const double zero = 0.;
        static const int lone = 1;

        int m  = num_rows(say);
        int n  = num_cols(sax);
        int ri = m-i;
        int rj = n-i-1;
        int ari = AMBIENT_IB-i-1;

        T* sayd = s_updated(say); int ldsay = say.num_rows();
        T* saxd = current(sax);   int ldsax = sax.num_rows();
        T* syd  = s_updated(sy);  int ldsy = sy.num_rows();
        T* sxd  = current(sx);    int ldsx = sx.num_rows();
        
        dgemv_("T", &ri, &ari, &one, &sayd[i + (i+1)*ldsay], &ldsay, &sayd[i+i*ldsay], &lone, &zero, &syd[i+1 + i*ldsy], &lone); // part of big gemv

        dgemv_("T", &ri, &i, &one, &sayd[i], &ldsay, &sayd[i+i*ldsay], &lone, &zero, &syd[i*ldsy], &lone);
        dgemv_("N", &rj, &i, &mone, &syd[ i+1 ], &ldsy, &syd[i*ldsy], &lone, &one, &syd[i+1 + i*ldsy], &lone);

        dgemv_("T", &ri, &i, &one, &sxd[i], &ldsx, &sayd[i + i*ldsay], &lone, &zero, &syd[ i*ldsy ], &lone);
        dgemv_("T", &i, &rj, &mone, &saxd[ (i+1)*ldsax ], &ldsax, &syd[i*ldsy], &lone, &one, &syd[ i+1 + i*ldsy], &lone);
    }

    template<typename T>
    void labrd_update_row<T>::c(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, matrix<T>& tp, matrix<T>& e, const int& i){
        static const double mone = -1.;
        static const double one = 1.;
        static const double zero = 0.;
        static const int lone = 1;

        int m   = num_rows(say);
        int n   = num_cols(sax);
        int ri  = m-i;
        int rj  = n-i-1;
        int rij = m-i-1;
        int r3  = i+1;

        T* sayd = current(say);     int ldsay = say.num_rows();
        T* saxd = s_updated(sax);   int ldsax = sax.num_rows();
        T* syd  = current(sy);      int ldsy = sy.num_rows();
        T* sxd  = s_updated(sx);    int ldsx = sx.num_rows();
        T* tpd  = s_updated(tp);
        T* ed   = s_updated(e);
        
        ambient::memptf<T, ambient::memcpy>(saxd, ldsax, dim2(i, i), 
                                            sayd, ldsay, dim2(i, i), 
                                            dim2( 1, ldsax-i ));
        
        dgemv_("T", &i, &rj, &mone, &saxd[(i+1)*ldsax], &ldsax, &sxd[i], &ldsx, &one, &saxd[ i + (i+1)*ldsax], &ldsax);
        dgemv_("N", &rj, &r3, &mone, &syd[ i+1 ], &ldsy, &saxd[i], &ldsax, &one, &saxd[i + (i+1)*ldsax], &ldsax);

        dlarfg_(&rj, &saxd[i + (i+1)*ldsax], &saxd[i + std::min(i+2, n-1)*ldsax], &ldsax, &tpd[i] );
        ed[i] = saxd[i + (i+1)*ldsax];
        saxd[i + (i+1)*ldsax] = 1.000;
    }

    template<typename T>
    void labrd_reduce_row<T>::c(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, const int& i){
        static const double mone = -1.;
        static const double one = 1.;
        static const double zero = 0.;
        static const int lone = 1;

        int m   = num_rows(say);
        int n   = num_cols(sax);
        int ri  = m-i;
        int rj  = n-i-1;
        int rij = m-i-1;
        int r3  = i+1;
        int ari = AMBIENT_IB-i-1;

        T* sayd = current(say);     int ldsay = say.num_rows();
        T* saxd = s_updated(sax);   int ldsax = sax.num_rows();
        T* syd  = current(sy);      int ldsy = sy.num_rows();
        T* sxd  = s_updated(sx);    int ldsx = sx.num_rows();
        
        dgemv_("T", &rj, &r3, &one, &syd[i+1], &ldsy, &saxd[i+(i+1)*ldsax], &ldsax, &zero, &sxd[i*ldsx], &lone);
        dgemv_("N", &rij, &r3, &mone, &sayd[i+1], &ldsay, &sxd[i*ldsx], &lone, &zero, &sxd[i+1+i*ldsx], &lone);

        dgemv_("N", &i, &rj, &one, &saxd[(i+1)*ldsax], &ldsax, &saxd[ i +(i+1)*ldsax], &ldsax, &zero, &sxd[i*ldsx], &lone);
        dgemv_("N", &rij, &i, &mone, &sxd[i+1], &ldsx, &sxd[i*ldsx], &lone, &one, &sxd[i+1+i*ldsx], &lone);

        dgemv_("N", &ari, &rj, &one, &saxd[i+1 + (i+1)*ldsax], &ldsax, &saxd[ i +(i+1)*ldsax], &ldsax, &one, &sxd[i+1 + i*ldsx], &lone); // part of big gemv
    }

    template<typename T, PLASMA_enum TR>
    void larfg<T,TR>::c(matrix<T>& a, matrix<T>& t, matrix<T>& d, const size_t& k){
        int lda;
        int n;
        T* alfa;
        T* x;
        
        T* ad = (T*)s_updated(a);
        T* td = (T*)s_updated(t);
        T* dd = (T*)s_updated(d);

        if(TR == PlasmaNoTrans){
            alfa = &ad[k + k*a.num_rows()];
            x = &ad[std::min(k+1, a.num_rows()-1)+k*a.num_rows()];
            n = a.num_rows()-k;
            lda = 1;
        }else{
            alfa = &ad[k + (k+1)*a.num_rows()];
            x = &ad[k + std::min(k+2, a.num_cols()-1)*a.num_rows()];
            n = a.num_cols()-k-1;
            lda = a.num_rows();
        }

        dlarfg_(&n, alfa, x, &lda, &td[k]);
        
        dd[k] = *alfa;
        *alfa = 1.00;
    }

    template<typename T>
    void gebd2<T>::c(matrix<T>& a, matrix<T>& d, matrix<T>& e, matrix<T>& tq, matrix<T>& tp){
        int m = a.num_rows();
        int n = a.num_cols();
        int lda = a.num_rows();
        int info;

        T* work = (T*)malloc(std::max(m,n)*sizeof(T));
        dgebd2_(&m, &n, (T*)s_updated(a), &lda, (T*)updated(d), (T*)updated(e), (T*)updated(tq), (T*)updated(tp), work, &info);
        free(work);
    }

    template<typename T>
    void gebrd<T>::c(matrix<T>& a, weak_view<T>& d, weak_view<T>& e, weak_view<T>& q, weak_view<T>& p){
        static int zero = 0;
        static int one  = 1;
        int m = q.num_rows();
        int n = p.num_rows();
        int lda = a.num_rows();
        int k = std::min(m,n);
        int lwork = -1;
        int info;
        
        T* ad = (T*)current(a);
        T* dd = (T*)updated(d);
        T* ed = (T*)updated(e);
        T* qd = (T*)updated(q);
        T* pd = (T*)updated(p);

        T* work = (T*)malloc(sizeof(T));
        T* tauq = (T*)malloc(sizeof(T)*k);
        T* taup = (T*)malloc(sizeof(T)*k);

        dgebrd_(&m, &n, ad, &lda, dd, ed, tauq, taup, work, &lwork, &info);
        lwork = (int)work[0];
        work = (T*)realloc(work, sizeof(T)*lwork);
        dgebrd_(&m, &n, ad, &lda, dd, ed, tauq, taup, work, &lwork, &info);

        T* ac = (T*)malloc(m*n*sizeof(T));
        std::memcpy(ac, ad, m*n*sizeof(T));

        lwork = -1;
        dorgbr_("Q",&m,&k,&n, ad, &lda, tauq, work, &lwork, &info);
        lwork = (int)work[0];
        work = (T*)realloc(work, sizeof(T)*lwork);
        dorgbr_("Q",&m,&k,&n, ad, &lda, tauq, work, &lwork, &info);


        lwork = -1;
        dorgbr_("P",&k,&n,&m, ac, &lda, taup, work, &lwork, &info);
        lwork = (int)work[0];
        work = (T*)realloc(work, sizeof(T)*lwork);
        dorgbr_("P",&k,&n,&m, ac, &lda, taup, work, &lwork, &info);

        for(int j = 0; j < n; ++j){
            std::memcpy(&pd[n*j], &ac[lda*j], sizeof(T)*k);
        }

        free(work); free(ac);
    }

    template<typename T>
    void gbbrd<T>::c(matrix<T>& a, weak_view<T>& d, weak_view<T>& e, weak_view<T>& q, weak_view<T>& p){
        static int zero = 0;
        static int one  = 1;
        int m = q.num_rows();
        int n = p.num_rows();
        int k = a.num_rows();
        int kl = (m <  n) ? k-1 : 0;
        int ku = (m >= n) ? k-1 : 0;
        int info;
        
        T* work = (T*)malloc(std::max(m,n)*2*sizeof(T));
        T* ad = (T*)current(a);
        T* dd = (T*)updated(d);
        T* ed = (T*)updated(e);
        T* qd = (T*)updated(q);
        T* pd = (T*)updated(p);

        dgbbrd_("B", &m, &n, &zero, &kl, &ku, ad, &k, dd, ed, 
                qd, &m, pd, &n, NULL, &one, work, &info);

        free(work);
    }

    template<typename T>
    void bdsqr<T>::c(matrix<T>& d, matrix<T>& e, matrix<T>& u, matrix<T>& v){
        static int zero = 0;
        static int one  = 1;
        int n = d.num_rows();
        int nv = v.num_cols();
        int mv = v.num_rows();
        int mu = u.num_rows();
        int info;
        
        T* work = (T*)malloc(n*4*sizeof(T));
        dbdsqr_("U", &n, &nv, &mu, &zero, (T*)s_updated(d), (T*)s_updated(e), 
                (T*)s_updated(v), &mv, (T*)s_updated(u), &mu, NULL, &one, work, &info); free(work);
        // Uncomment for dc numerically loose algorithm:
        // LAPACKE_dbdsdc(102, 'U', 'N', n, (T*)s_updated(d), (T*)s_updated(e),  (T*)s_updated(u), one, (T*)s_updated(v), one, NULL, NULL);
    }

    template<typename T, PLASMA_enum UL>
    void copy_band<T,UL>::c(const matrix<T>& src, matrix<T>& dst, const size_t& dj){
        T* sd = (T*)current(src);
        T* dd = (T*)s_updated(dst);
        size_t ldd = dst.num_rows();
        size_t m = src.num_rows();
        size_t n = src.num_cols();
        size_t offset = std::min(ldd-1,(size_t)AMBIENT_IB);

        dd += dj*ldd;
        if(UL == PlasmaUpper){
            for(int j = 0; j < n; ++j)
            for(int i = 0; i <= j && i < m; ++i)
            dd[j*ldd+i-j+offset] = sd[i+m*j]; 
        }else{
            for(int j = 0; j < n; ++j)
            for(int i = j; i < m; ++i)
            dd[j*ldd+i-j] = sd[i+m*j]; 
        }
    }

    template<int ADD, class VA, class VB, class VC, class VF>
    void gemv_scale<ADD, VA, VB, VC, VF>::c(const matrix<T>& a, const size_t& aoffset, 
                                            const matrix<T>& b, const size_t& boffset,
                                                  matrix<T>& c, const size_t& coffset,
                                            const matrix<T>& f, const size_t& foffset,
                                            const size_t& rows, const size_t& cols)
    {
        T* ad = current(a);
        T* bd = current(b);
        T* fd = current(f);
        T* cd = s_updated(c);
        int lda = num_rows(a);
        int ldb = VB::inc(b);
        int ldc = VC::inc(c);
        int m = rows;
        int n = cols;
        static const double one = 1.;
        const double* alfa = &fd[foffset];
        const double* beta = (ADD == 1) ? &one : alfa;

        if(*VA::code() == 'T') std::swap(m,n);
        dgemv_(VA::code(), &m, &n, alfa, &ad[aoffset], &lda, &bd[boffset], &ldb, beta, &cd[coffset], &ldc);
    }

    template<int alfa, int beta, class ViewA, class ViewB, class ViewC>
    void gemv<alfa, beta, ViewA, ViewB, ViewC>::c(const matrix<T>& a, const size_t& aoffset, 
                  const matrix<T>& b, const size_t& boffset,
                        matrix<T>& c, const size_t& coffset,
                  const size_t& rows, const size_t& cols)
    {
        T* ad = current(a);
        T* bd = current(b);
        T* cd = s_updated(c);
        int lda = num_rows(a);
        int ldb = ViewB::inc(b);
        int ldc = ViewC::inc(c);
        static const double salfa(alfa); 
        static const double sbeta(beta);
        int m = rows;
        int n = cols;

        if(*ViewA::code() == 'T') std::swap(m,n);
        dgemv_(ViewA::code(), &m, &n, &salfa, &ad[aoffset], &lda, &bd[boffset], &ldb, &sbeta, &cd[coffset], &ldc);
    }

} } }
