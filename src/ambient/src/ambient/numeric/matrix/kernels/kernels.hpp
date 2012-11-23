#ifndef AMBIENT_NUMERIC_MATRIX_KERNELS
#define AMBIENT_NUMERIC_MATRIX_KERNELS
#define pin this->pin
#define AMBIENT_IB 512
#define PLASMA_IB 64

namespace ambient { namespace numeric { namespace kernels {

    using ambient::numeric::matrix;
    using ambient::numeric::weak_view;

    // kernels required by QR //

    template<typename T>
    struct geqrt : public kernel< geqrt<T> > 
    {
        typedef void(geqrt::*F)(matrix<T>&, matrix<T>&);

        inline void l(matrix<T>& a, matrix<T>& t){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(t));
        }

        inline void c(matrix<double>& a, matrix<double>& t){
            __A_TIME_C("ambient_geqrt_c_kernel"); 
            T* tau  = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB>(); 
            T* work = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
            CORE_dgeqrt(a.num_rows(), a.num_cols(), PLASMA_IB,
                        (T*)s_updated(a), a.num_rows(),
                        (T*)updated(t),   t.num_rows(), 
                        tau, work);
            __A_TIME_C_STOP
        }

        inline void c(matrix<std::complex<double> >& a, matrix<std::complex<double> >& t){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    template<typename T, PLASMA_enum TR>
    struct ormqr : public kernel< ormqr<T,TR> > 
    {
        typedef void(ormqr::*F)(const size_t&, const matrix<T>&, const matrix<T>&, matrix<T>&);

        inline void l(const size_t& k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(t));
            assign(current(c));
        }

        inline void c(const size_t& k, const matrix<double>& a, const matrix<double>& t, matrix<double>& c){
            __A_TIME_C("ambient_ormqr_c_kernel"); 
            T* work = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
            CORE_dormqr(PlasmaLeft, TR,
                        c.num_rows(), c.num_cols(), k, PLASMA_IB,
                        (T*)current(a), a.num_rows(),
                        (T*)current(t), t.num_rows(),
                        (T*)s_updated(c), c.num_rows(),
                        work, AMBIENT_IB);
            __A_TIME_C_STOP
        }

        inline void c(const size_t& k, const matrix<std::complex<double> >& a, const matrix<std::complex<double> >& t, matrix<std::complex<double> >& c){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    template<typename T>
    struct tsqrt : public kernel< tsqrt<T> > 
    {
        typedef void(tsqrt::*F)(matrix<T>&, matrix<T>&, matrix<T>&);

        inline void l(matrix<T>& a1, matrix<T>& a2, matrix<T>& t){
            pin(current(a1)); //if(!ctxt.involved()) return;
            assign(current(a2));
            assign(current(t));
        }

        inline void c(matrix<double>& a1, matrix<double>& a2, matrix<double>& t){
            __A_TIME_C("ambient_tsqrt_c_kernel"); 
            T* tau  = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB>();
            T* work = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
            CORE_dtsqrt(a2.num_rows(), a2.num_cols(), PLASMA_IB,
                        (T*)s_updated(a1), a1.num_rows(),
                        (T*)s_updated(a2), a2.num_rows(),
                        (T*)updated(t),     t.num_rows(),
                        tau, work);
            __A_TIME_C_STOP
        }

        inline void c(matrix<std::complex<double> >& a1, matrix<std::complex<double> >& a2, matrix<std::complex<double> >& t){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    template<typename T, PLASMA_enum TR>
    struct tsmqr : public kernel< tsmqr<T,TR> > 
    {
        typedef void(tsmqr::*F)(const size_t&, matrix<T>&, matrix<T>&, const matrix<T>&, const matrix<T>&);

        inline void l(const size_t& k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t){
            pin(current(a1)); //if(!ctxt.involved()) return;
            assign(current(a2));
            assign(current(v));
            assign(current(t));
        }

        inline void c(const size_t& k, matrix<double>& a1, matrix<double>& a2, const matrix<double>& v, const matrix<double>& t){
            __A_TIME_C("ambient_tsmqr_c_kernel"); 
            T* work = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
            CORE_dtsmqr(PlasmaLeft, TR,
                        AMBIENT_IB, a1.num_cols(), a2.num_rows(), a2.num_cols(), k, PLASMA_IB,
                        (T*)s_updated(a1), a1.num_rows(),
                        (T*)s_updated(a2), a2.num_rows(),
                        (T*)current(v), v.num_rows(),
                        (T*)current(t), t.num_rows(),
                        work, PLASMA_IB);
            __A_TIME_C_STOP
        }

        inline void c(const size_t& k, matrix<std::complex<double> >& a1, matrix<std::complex<double> >& a2, const matrix<std::complex<double> >& v, const matrix<std::complex<double> >& t){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    // kernels required by LQ //

    template<typename T>
    struct gelqt : public kernel< gelqt<T> > 
    {
        typedef void(gelqt::*F)(matrix<T>&, matrix<T>&);

        inline void l(matrix<T>& a, matrix<T>& t){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(t));
        }

        inline void c(matrix<double>& a, matrix<double>& t){
            __A_TIME_C("ambient_gelqt_c_kernel"); 
            T* tau  = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB>();
            T* work = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
            CORE_dgelqt(a.num_rows(), a.num_cols(), PLASMA_IB,
                        (T*)s_updated(a), a.num_rows(), 
                        (T*)updated(t),   t.num_rows(),
                        tau, work);
            __A_TIME_C_STOP
        }

        inline void c(matrix<std::complex<double> >& a, matrix<std::complex<double> >& t){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    template<typename T, PLASMA_enum TR>
    struct ormlq : public kernel< ormlq<T,TR> > 
    {
        typedef void(ormlq::*F)(const size_t&, const matrix<T>&, const matrix<T>&, matrix<T>&);

        inline void l(const size_t& k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(t));
            assign(current(c));
        }

        inline void c(const size_t& k, const matrix<double>& a, const matrix<double>& t, matrix<double>& c){
            __A_TIME_C("ambient_ormlq_c_kernel"); 
            T* work = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
            CORE_dormlq(PlasmaRight, TR,
                        c.num_rows(), c.num_cols(), k, PLASMA_IB,
                        (T*)current(a), a.num_rows(),
                        (T*)current(t), t.num_rows(),
                        (T*)s_updated(c), c.num_rows(),
                        work, AMBIENT_IB);
            __A_TIME_C_STOP
        }

        inline void c(const size_t& k, const matrix<std::complex<double> >& a, const matrix<std::complex<double> >& t, matrix<std::complex<double> >& c){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    template<typename T>
    struct tslqt : public kernel< tslqt<T> > 
    {
        typedef void(tslqt::*F)(matrix<T>&, matrix<T>&, matrix<T>&);

        inline void l(matrix<T>& a1, matrix<T>& a2, matrix<T>& t){
            pin(current(a1)); //if(!ctxt.involved()) return;
            assign(current(a2));
            assign(current(t));
        }

        inline void c(matrix<double>& a1, matrix<double>& a2, matrix<double>& t){
            __A_TIME_C("ambient_tslqt_c_kernel"); 
            T* tau  = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB>();
            T* work = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
            CORE_dtslqt(a2.num_rows(), a2.num_cols(), PLASMA_IB,
                        (T*)s_updated(a1), a1.num_rows(),
                        (T*)s_updated(a2), a2.num_rows(),
                        (T*)updated(t),     t.num_rows(),
                        tau, work);
            __A_TIME_C_STOP
        }

        inline void c(matrix<std::complex<double> >& a1, matrix<std::complex<double> >& a2, matrix<std::complex<double> >& t){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    template<typename T, PLASMA_enum TR>
    struct tsmlq : public kernel< tsmlq<T,TR> > 
    {
        typedef void(tsmlq::*F)(const size_t&, matrix<T>&, matrix<T>&, const matrix<T>&, const matrix<T>&);

        inline void l(const size_t& k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t){
            pin(current(a1)); //if(!ctxt.involved()) return;
            assign(current(a2));
            assign(current(v));
            assign(current(t));
        }

        inline void c(const size_t& k, matrix<double>& a1, matrix<double>& a2, const matrix<double>& v, const matrix<double>& t){
            __A_TIME_C("ambient_tsmlq_c_kernel"); 
            T* work = (T*)ambient::bulk_pool.get<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
            CORE_dtsmlq(PlasmaRight, TR,
                        a1.num_rows(), AMBIENT_IB, a2.num_rows(), a2.num_cols(), k, PLASMA_IB,
                        (T*)s_updated(a1), a1.num_rows(),
                        (T*)s_updated(a2), a2.num_rows(),
                        (T*)current(v), v.num_rows(),
                        (T*)current(t), t.num_rows(),
                        work, AMBIENT_IB);
            __A_TIME_C_STOP
        }

        inline void c(const size_t& k, matrix<std::complex<double> >& a1, matrix<std::complex<double> >& a2, const matrix<std::complex<double> >& v, const matrix<std::complex<double> >& t){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    template<typename T, PLASMA_enum UL, size_t OFF>
    struct laset2 : public kernel< laset2<T,UL,OFF> > 
    {
        typedef void(laset2::*F)(matrix<T>&, const T&);

        inline void l(matrix<T>& a, const T&){
            pin(current(a)); //if(!ctxt.involved()) return;
        }

        inline void c(matrix<double>& a, const double& alfa){
            __A_TIME_C("ambient_laset2_c_kernel"); 
            T* ad = (T*)s_updated(a);
            CORE_dlaset2(UL, a.num_rows()-OFF, a.num_cols()-OFF, alfa, ad + OFF*a.num_rows(), a.num_rows());
            __A_TIME_C_STOP
        }

        inline void c(matrix<std::complex<double> >& a, const std::complex<double>& alfa){
            printf("NOT IMPLEMENTED!\n");
        }
    };

    template<int ADD, class VA, class VB, class VC, class VF>
    struct gemv_scale : public kernel< gemv_scale<ADD, VA, VB, VC, VF> > 
    { // gs
        typedef typename VC::value_type T;
        typedef void(gemv_scale::*F)(const matrix<T>&, const size_t&, 
                                     const matrix<T>&, const size_t&, 
                                           matrix<T>&, const size_t&, 
                                     const matrix<T>&, const size_t&, 
                                     const size_t&, const size_t&);

        inline void l(const matrix<T>& a, const size_t& aoffset, 
                      const matrix<T>& b, const size_t& boffset,
                            matrix<T>& c, const size_t& coffset,
                      const matrix<T>& f, const size_t& foffset,
                      const size_t& rows, const size_t& cols){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b));
            assign(current(c));
            assign(current(f));
        }

        inline void c(const matrix<T>& a, const size_t& aoffset, 
                      const matrix<T>& b, const size_t& boffset,
                            matrix<T>& c, const size_t& coffset,
                      const matrix<T>& f, const size_t& foffset,
                      const size_t& rows, const size_t& cols){
            __A_TIME_C("ambient_gemv_scale_c_kernel"); 
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
            __A_TIME_C_STOP
        }
    };

    template<int alfa, int beta, class ViewA, class ViewB, class ViewC>
    struct gemv : public kernel< gemv<alfa, beta, ViewA, ViewB, ViewC> > 
    { // gs
        typedef typename ViewC::value_type T;
        typedef void(gemv::*F)(const matrix<T>&, const size_t&, 
                               const matrix<T>&, const size_t&, 
                                     matrix<T>&, const size_t&, 
                               const size_t&, const size_t&);

        inline void l(const matrix<T>& a, const size_t& aoffset, 
                      const matrix<T>& b, const size_t& boffset,
                            matrix<T>& c, const size_t& coffset,
                      const size_t& rows, const size_t& cols){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b));
            assign(current(c));
        }

        inline void c(const matrix<T>& a, const size_t& aoffset, 
                      const matrix<T>& b, const size_t& boffset,
                            matrix<T>& c, const size_t& coffset,
                      const size_t& rows, const size_t& cols){
            __A_TIME_C("ambient_gemv_c_kernel"); 
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
            __A_TIME_C_STOP
        }
    };

    template<class ViewA, class ViewB, typename T>
    struct gemm : public kernel< gemm<ViewA, ViewB, T> > 
    { // gs
        typedef void(gemm::*F)(const matrix<T>&, const matrix<T>&, weak_view<T>&);

        inline void l(const matrix<T>& a, const matrix<T>& b, weak_view<T>& c){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b));
            assign(current(c));
        }

        inline void c(const matrix<double>& a, const matrix<double>& b, weak_view<double>& c){
            __A_TIME_C("ambient_gemm_general_c_kernel"); 
            if(!current(a).valid() || !current(b).valid()){
                (T*)p_updated(c);
                return;
            }

            double* ad = current(a);
            double* bd = current(b);
            double* cd = updated(c);
            int m = ViewA::rows(a);
            int k = ViewA::cols(a);
            int n = ViewB::cols(b);
            int lda = __a_get_dim(a).y;
            int ldb = __a_get_dim(b).y;
            int ldc = __a_get_dim(c).y;
            static const double alpha(1.0); 
            static const double beta(0.0);
            dgemm_(ViewA::code(), ViewB::code(), &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
            __A_TIME_C_STOP
        }
        inline void c(const matrix<std::complex<double> >& a, const matrix<std::complex<double> >& b, weak_view<std::complex<double> >& c){
            __A_TIME_C("ambient_gemm_general_c_kernel"); 
            T* ad   = current(a);
            T* bd   = current(b);
            T* cd   = updated(c);
            int m   = __a_get_dim(a).y;
            int n   = __a_get_dim(b).x;
            int k   = __a_get_dim(b).y;
            T alpha(1.0); 
            T beta(0.0);
            __a_gemm("N","N", &m, &n, &k, &alpha, ad, &m, bd, &k, &beta, cd, &m);
            __A_TIME_C_STOP
        }
    };
        
    template<class ViewB, typename T, typename D>
    struct gemm_diagonal_lhs : public kernel< gemm_diagonal_lhs<ViewB,T,D> > 
    {
        typedef void (gemm_diagonal_lhs::*F)(const matrix<D>&, const matrix<T>&, weak_view<T>&);

        inline void l(const matrix<D>& a_diag, const matrix<T>& b, weak_view<T>& c){
            pin(current(b)); //if(!ctxt.involved()) return;
            assign(current(a_diag));
            assign(current(c));
        }

        inline void c(const matrix<D>& a_diag, const matrix<T>& b, weak_view<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_diagonal_lhs_c_kernel"); 
            int sizey = __a_get_dim(a_diag).y;
            int size = __a_get_dim(b).x;
            int ONE  = 1;
            D* bd = current(b);
            D* cd = p_updated(c);
            D* alpha = current(a_diag);
        
            for(int k = 0 ; k < sizey; k++){
        	     axpy(&size, &alpha[k], &bd[k], &sizey, &cd[k], &sizey);
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T, typename D>
    struct gemm_diagonal_lhs<transpose_view<matrix<T> >,T,D> : public kernel< gemm_diagonal_lhs<transpose_view<matrix<T> >,T,D> > 
    {
        typedef void (gemm_diagonal_lhs::*F)(const matrix<D>&, const matrix<T>&, weak_view<T>&);

        inline void l(const matrix<D>& a_diag, const matrix<T>& b, weak_view<T>& c){
            pin(current(b)); //if(!ctxt.involved()) return;
            assign(current(a_diag));
            assign(current(c));
        }

        inline void c(const matrix<D>& a_diag, const matrix<T>& b, weak_view<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_diagonal_lhs_c_kernel"); 
            printf("Special DIAGONAL!\n");
            size_t sizex = __a_get_dim(b).x;
            int size  = __a_get_dim(a_diag).y;
            int ONE  = 1;
            D* bd = current(b);
            D* cd = p_updated(c);
            D* alpha = current(a_diag);
        
            for(int k = 0 ; k < sizex; k++){
        	     axpy(&size, &alpha[k], &bd[k*size], &ONE, &cd[k], &size);
            }
            __A_TIME_C_STOP
        }
    };
        
    template<class ViewA, typename T, typename D>
    struct gemm_diagonal_rhs : public kernel< gemm_diagonal_rhs<ViewA,T,D> > 
    {
        typedef void (gemm_diagonal_rhs::*F)(const matrix<T>&, const matrix<D>&, weak_view<T>&);

        inline void l(const matrix<T>& a, const matrix<D>& b_diag, weak_view<T>& c){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b_diag));
            assign(current(c));
        }

        inline void c(const matrix<T>& a, const matrix<D>& b_diag, weak_view<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_diagonal_rhs_c_kernel"); 
            size_t sizex = __a_get_dim(b_diag).y;
            int size = __a_get_dim(a).y; // for the case of complex
            int ONE = 1;
            D* ad = current(a);
            D* cd = p_updated(c);
        	D* alpha = current(b_diag);
        
            for(int k = 0 ; k < sizex; k++){
        	    axpy(&size, &alpha[k], &ad[k*size], &ONE, &cd[k*size], &ONE);
            }
            __A_TIME_C_STOP
        }
    };

    template<typename T, typename D>
    struct gemm_diagonal_rhs<transpose_view<matrix<T> >,T,D> : public kernel< gemm_diagonal_rhs<transpose_view<matrix<T> >,T,D> > 
    {
        typedef void (gemm_diagonal_rhs::*F)(const matrix<T>&, const matrix<D>&, weak_view<T>&);

        inline void l(const matrix<T>& a, const matrix<D>& b_diag, weak_view<T>& c){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b_diag));
            assign(current(c));
        }

        inline void c(const matrix<T>& a, const matrix<D>& b_diag, weak_view<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_diagonal_rhs_c_kernel"); 
            printf("Special DIAGONAL!\n");
            int sizey = __a_get_dim(b_diag).y;
            int size = __a_get_dim(a).x;
            int ONE = 1;
            D* ad = current(a);
            D* cd = p_updated(c);
        	D* alpha = current(b_diag);
        
            for(int k = 0 ; k < sizey; k++){
        	    axpy(&size, &alpha[k], &ad[k], &sizey, &cd[k*size], &ONE);
            }
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct copy : public kernel< copy<T> > 
    { // gs
        typedef void(copy::*F)(weak_view<T>&, const matrix<T>&);

        inline void l(weak_view<T>& ac, const matrix<T>& a){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(ac));
        }

        inline void c(weak_view<T>& ac, const matrix<T>& a){
            __A_TIME_C("ambient_copy_c_kernel"); 
            if(!current(a).valid()) (T*)p_updated(ac);
            else __a_refresh<T>(updated(ac), current(a), __a_sizeof(a));
            __A_TIME_C_STOP
        }
    };

    template<typename T, PLASMA_enum UL>
    struct copy_band : public kernel< copy_band<T,UL> > 
    { // gs
        typedef void(copy_band::*F)(matrix<T>&, const size_t&, const matrix<T>&); 

        inline void l(matrix<T>& dst, const size_t& dj, const matrix<T>& src)
        {
            pin(current(src)); //if(!ctxt.involved()) return;
            assign(current(dst));
        }

        inline void c(matrix<T>& dst, const size_t& dj, const matrix<T>& src)
        {
            __A_TIME_C("ambient_copy_band_c_kernel"); 
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
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct copy_rt : public kernel< copy_rt<T> > 
    { // gs
        typedef void(copy_rt::*F)(weak_view<T>&, const matrix<T>&);

        inline void l(weak_view<T>& t, const matrix<T>& a){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(t));
        }

        inline void c(weak_view<T>& t, const matrix<T>& a){
            __A_TIME_C("ambient_copy_rt_kernel"); 
            T* ad  = current(a);
            T* td  = p_updated(t);
            size_t sda = a.num_cols();
            size_t lda = a.num_rows();
            size_t ldt = t.num_rows();

            for(int j = 0; j < sda; ++j)
            for(int i = 0; i <= j && i < ldt; ++i)
            td[i+ldt*j] = ad[i+lda*j]; 
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct copy_lt : public kernel< copy_lt<T> > 
    { // gs
        typedef void(copy_lt::*F)(weak_view<T>&, const matrix<T>&);

        inline void l(weak_view<T>& t, const matrix<T>& a){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(t));
        }

        inline void c(weak_view<T>& t, const matrix<T>& a){
            __A_TIME_C("ambient_copy_lt_kernel"); 
            T* ad  = current(a);
            T* td  = p_updated(t);
            size_t sdt = t.num_cols();
            size_t lda = a.num_rows();
            size_t ldt = t.num_rows();

            for(int j = 0; j < sdt; ++j)
            for(int i = j; i < lda; ++i)
            td[i+ldt*j] = ad[i+lda*j]; 
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct copy_partial : public kernel< copy_partial<T> > 
    { // gs
        typedef void(copy_partial::*F)(matrix<T>&, const size_t&, const size_t&, 
                                       const matrix<T>&, const size_t&, const size_t&, 
                                       const size_t&, const size_t&);
                                           

        inline void l(matrix<T>& dst, const size_t& di, const size_t& dj, 
                      const matrix<T>& src, const size_t& si, const size_t& sj,
                      const size_t& m, const size_t& n)
        {
            pin(current(src)); //if(!ctxt.involved()) return;
            assign(current(dst));
        }

        inline void c(matrix<T>& dst, const size_t& di, const size_t& dj, 
                      const matrix<T>& src, const size_t& si, const size_t& sj,
                      const size_t& m, const size_t& n)
        {
            __A_TIME_C("ambient_copy_partial_c_kernel"); 
            T* sd = current(src);
            T* dd = m*n < __a_get_dim(dst).square() ? (T*)s_updated(dst) : (T*)updated(dst);
            __a_memptf_r<T, __a_memcpy>(dd, __a_get_dim(dst).y, dim2(dj, di), 
                                        sd, __a_get_dim(src).y, dim2(sj, si), 
                                        dim2( n, m ));
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct copy_s : public kernel< copy_s<T> > 
    { // gs
        typedef void(copy_s::*F)(matrix<T>&, const size_t&, const size_t&, 
                                 const matrix<T>&, const size_t&, const size_t&, 
                                 const matrix<T>&, const size_t&, const size_t&, 
                                 const size_t&, const size_t&);
                                           
        inline void l(matrix<T>& dst, const size_t& di, const size_t& dj, 
                      const matrix<T>& src, const size_t& si, const size_t& sj,
                      const matrix<T>& alfa, const size_t& ai, const size_t& aj,
                      const size_t& m, const size_t& n)
        {
            pin(current(src)); //if(!ctxt.involved()) return;
            assign(current(dst));
            assign(current(alfa));
        }

        inline void c(matrix<T>& dst, const size_t& di, const size_t& dj, 
                      const matrix<T>& src, const size_t& si, const size_t& sj,
                      const matrix<T>& alfa, const size_t& ai, const size_t& aj,
                      const size_t& m, const size_t& n)
        {
            __A_TIME_C("ambient_copy_s_c_kernel"); 
            T factor = ((T*)current(alfa))[ai + aj*__a_get_dim(alfa).y];
            T* sd = current(src);
            T* dd = m*n < __a_get_dim(dst).square() ? (T*)s_updated(dst) : (T*)updated(dst);
            __a_memptf_r<T, __a_memscal>(dd, __a_get_dim(dst).y, dim2(dj, di), 
                                         sd, __a_get_dim(src).y, dim2(sj, si), 
                                         dim2( n, m ), factor);
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct copy_sa : public kernel< copy_sa<T> > 
    { // gs
        typedef void(copy_sa::*F)(matrix<T>&, const size_t&, const size_t&, 
                                  const matrix<T>&, const size_t&, const size_t&, 
                                  const matrix<T>&, const size_t&, const size_t&, 
                                  const size_t&, const size_t&);
                                           
        inline void l(matrix<T>& dst, const size_t& di, const size_t& dj, 
                      const matrix<T>& src, const size_t& si, const size_t& sj,
                      const matrix<T>& alfa, const size_t& ai, const size_t& aj,
                      const size_t& m, const size_t& n)
        {
            pin(current(src)); //if(!ctxt.involved()) return;
            assign(current(dst));
            assign(current(alfa));
        }

        inline void c(matrix<T>& dst, const size_t& di, const size_t& dj, 
                      const matrix<T>& src, const size_t& si, const size_t& sj,
                      const matrix<T>& alfa, const size_t& ai, const size_t& aj,
                      const size_t& m, const size_t& n)
        {
            __A_TIME_C("ambient_copy_sa_c_kernel"); 
            T factor = ((T*)current(alfa))[ai + aj*__a_get_dim(alfa).y];
            __a_memptf_r<T, __a_memscala>(s_updated(dst), __a_get_dim(dst).y, dim2(dj, di), 
                                          current(src), __a_get_dim(src).y, dim2(sj, si), 
                                          dim2( n, m ), factor);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct trace : public kernel< trace<T> > 
    {
        typedef void (trace::*F)(const matrix<T>&, future<T>&);

        inline void l(const matrix<T>& a, future<T>& trace){
            pin(current(a));
        }

        inline void c(const matrix<T>& a, future<T>& trace){
            // gs
            __A_TIME_C("ambient_trace_c_kernel"); 
            size_t m = __a_get_dim(a).y;
            size_t n = __a_get_dim(a).x;
            T* ad = current(a);
        
            size_t sizex = std::min(n,m);
            for(size_t jj = 0; jj < sizex; jj++){
                trace.get_value() += ad[jj + jj*m];
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scalar_norm : public kernel< scalar_norm<T> > 
    {// gs
        typedef void (scalar_norm::*F)(const matrix<T>&, future<double>&);

        inline void l(const matrix<T>& a, future<double>& norm){
            pin(current(a));
        }

        inline void c(const matrix<T>& a, future<double>& norm){
            __A_TIME_C("ambient_scalar_norm_c_kernel"); 
            T* ad = current(a);
            norm.get_value() = alps::numeric::real(__a_dot(ad, ad, __a_get_dim(a).square()));
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct overlap : public kernel< overlap<T> > 
    { // gs
        typedef void (overlap::*F)(const matrix<T>&, const matrix<T>&, future<T>&);

        inline void l(const matrix<T>& a, const matrix<T>& b, future<T>& overlap){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b));
        }

        inline void c(const matrix<T>& a, const matrix<T>& b, future<T>& overlap){
            __A_TIME_C("ambient_scalar_overlap_c_kernel"); 
            T* ad = current(a);
            T* bd = current(b);
            overlap.get_value() = __a_dot(ad, bd, __a_get_dim(a).square());
            __A_TIME_C_STOP
        }
    };

        
    template<int alfa, typename T>
    struct add_vectors : public kernel< add_vectors<alfa, T> > 
    { // gs
        typedef void (add_vectors::*F)(matrix<T>&, const size_t&, const matrix<T>&, const size_t&, const size_t&);

        inline void l(matrix<T>& a, const size_t& aoffset, const matrix<T>& b, const size_t& boffset, const size_t& size){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b));
        }

        inline void c(matrix<T>& a, const size_t& aoffset, const matrix<T>& b, const size_t& boffset, const size_t& size){
            __A_TIME_C("ambient_add_vectors_c_kernel");
            T* ad = &((T*)current(a))[aoffset];
            T* bd = &((T*)current(b))[boffset];
            T* ar = &((T*)updated(a))[aoffset];

            for(int k = 0; k < size; k++) 
                ar[k] = alfa*ad[k] + bd[k];
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct add : public kernel< add<T> > 
    { // gs
        typedef void (add::*F)(matrix<T>&, const matrix<T>&);

        inline void l(matrix<T>& a, const matrix<T>& b){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b));
        }

        inline void c(matrix<T>& a, const matrix<T>& b){
            __A_TIME_C("ambient_add_c_kernel");
            T* ad = current(a);
            T* bd = current(b);
            T* ar = updated(a);

            int size = __a_get_dim(a).square();
            for(int k = 0; k < size; k++) 
                ar[k] = ad[k] + bd[k];
            __A_TIME_C_STOP
        }
    };

        
    template<typename T>
    struct sub : public kernel< sub<T> > 
    { // gs
        typedef void (sub::*F)(matrix<T>&, const matrix<T>&);

        inline void l(matrix<T>& a, const matrix<T>& b){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(b));
        }

        inline void c(matrix<T>& a, const matrix<T>& b){
            __A_TIME_C("ambient_sub_c_kernel"); 
            T* ad = current(a);
            T* bd = current(b);
            T* ar = updated(a);

            int size = __a_get_dim(a).square();
            for(int k = 0; k < size; k++) 
                ar[k] = ad[k] - bd[k];
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scale : public kernel< scale<T> > 
    { // gs
        typedef void (scale::*F)(matrix<T>&, const future<T>&);

        inline void l(matrix<T>& a, const future<T>& t){
            pin(current(a));
        }

        inline void c(matrix<T>& a, const future<T>& t){
            __A_TIME_C("ambient_scale_c_kernel"); 
            T* ad = current(a);
            T* ar = updated(a);
            int size = __a_get_dim(a).square();
            for(int k=0; k < size; k++) 
                ar[k] = ad[k] * t.get_value();
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scale_offset : public kernel< scale_offset<T> > 
    { // gs
        typedef void (scale_offset::*F)(matrix<T>&, const size_t&, const size_t&, const matrix<T>&, const size_t&);

        inline void l(matrix<T>& a, const size_t& ai, const size_t& aj, const matrix<T>& alfa, const size_t& alfai){
            pin(current(a));
            assign(current(alfa));
        }

        inline void c(matrix<T>& a, const size_t& ai, const size_t& aj, const matrix<T>& alfa, const size_t& alfai){
            __A_TIME_C("ambient_scale_offset_c_kernel"); 
            int m = num_rows(a);
            T* ad = &((T*)s_updated(a))[aj*m];
            T factor = ((T*)current(alfa))[alfai];
            for(int k = ai; k < m; k++) ad[k] *= factor;
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scale_inverse : public kernel< scale_inverse<T> > 
    { // gs
        typedef void (scale_inverse::*F)(matrix<T>&, const future<T>&);

        inline void l(matrix<T>& a, const future<T>& t){
            pin(current(a));
        }

        inline void c(matrix<T>& a, const future<T>& t){
            __A_TIME_C("ambient_scale_inverse_c_kernel"); 
            T* ad = current(a);
            T* ar = updated(a);
            int size = __a_get_dim(a).square();
            for(int k=0; k < size; k++) 
                ar[k] = ad[k] / t.get_value();
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct sqrt_diagonal : public kernel< sqrt_diagonal<T> > 
    {
        typedef void (sqrt_diagonal::*F)(matrix<T>&);

        inline void l(matrix<T>& a){
            pin(current(a));
        }

        inline void c(matrix<T>& a){
            __A_TIME_C("ambient_sqrt_diagonal_c_kernel"); 
            size_t size = __a_get_dim(a).y;
            T* ad = current(a);
            T* ar = updated(a);
            for(size_t i = 0; i < size; ++i) ar[i] = std::sqrt(ad[i]);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct exp_diagonal : public kernel< exp_diagonal<T> > 
    {
        typedef void (exp_diagonal::*F)(matrix<T>&, const T&);

        inline void l(matrix<T>& a, const T& alfa){
            pin(current(a));
        }

        inline void c(matrix<T>& a, const T& alfa){
            __A_TIME_C("ambient_exp_diagonal_c_kernel"); 
            size_t size = __a_get_dim(a).y;
            T* ad = current(a);
            T* ar = updated(a);
            for(size_t i = 0; i < size; ++i) ar[i] = std::exp(alfa*ad[i]);
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct transpose_out : public kernel< transpose_out<T> > 
    { // gs
        typedef void (transpose_out::*F)(const matrix<T>&, weak_view<T>&);

        inline void l(const matrix<T>& a, weak_view<T>& t){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(t));
        }

        inline void c(const matrix<T>& a, weak_view<T>& t){
            __A_TIME_C("ambient_transpose_out_c_kernel"); 
            T* od = current(a);
            T* td = updated(t);
            int m = __a_get_dim(a).y;
            int n = __a_get_dim(a).x;

            for(int i = 0; i < m; i++){
                for(int j = 0; j < n; j++) *td++ = od[j*m];
                od++;
            }
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct resize : public kernel< resize<T> > 
    {
        typedef void (resize::*F)(weak_view<T>&, const matrix<T>&, const size_t&, const size_t&);

        inline void l(weak_view<T>& r, const matrix<T>& a, const size_t& m, const size_t& n){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(r));
        }

        inline void c(weak_view<T>& r, const matrix<T>& a, const size_t& m, const size_t& n){
            __A_TIME_C("ambient_resize_c_kernel");
            T* dd = m*n == __a_get_dim(r).square() ? (T*)updated(r) : (T*)p_updated(r);
            __a_memptf_r<T, __a_memcpy>(dd, __a_get_dim(r).y, dim2(0,0),
                                        current(a), __a_get_dim(a).y, dim2(0,0), dim2(n, m)); 
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct init_identity : public kernel< init_identity<T> > 
    {
        typedef void (init_identity::*F)(weak_view<T>&);

        inline void l(weak_view<T>& a){
            pin(current(a));
        }

        inline void c(weak_view<T>& a){
            __A_TIME_C("ambient_init_identity_c_kernel"); 
            size_t n = __a_get_dim(a).x;
            size_t m = __a_get_dim(a).y;
            T* ad = p_updated(a);

            size_t sizex = std::min(m,n); // respecting borders
            for(size_t jj = 0; jj < sizex; ++jj) ad[jj + m*jj] = 1.;
            __A_TIME_C_STOP
        }
    };
       
    template<typename T>
    struct init_random : public kernel< init_random<T> > 
    {
        typedef void (init_random::*F)(weak_view<T>&);
     
        template<typename T> inline void randomize(T& a){ 
            a = drand48();
        }
        template<typename T> inline void randomize(std::complex<T>& a){
            a.real(drand48());
            a.imag(drand48());
        }

        inline void l(weak_view<T>& a){
            pin(current(a));
        }
        
        inline void c(weak_view<T>& a){
            __A_TIME_C("ambient_init_random_c_kernel"); 
            size_t size = __a_get_dim(a).square();
            T* ad = updated(a);
            for(size_t i = 0; i < size; ++i) randomize(ad[i]);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct init_value : public kernel< init_value<T> > 
    {
        typedef void (init_value::*F)(weak_view<T>&, const T&);

        inline void l(weak_view<T>& a, const T& value){
            pin(current(a));
        }

        inline void c(weak_view<T>& a, const T& value){
            __A_TIME_C("ambient_init_value_c_kernel"); 
            size_t size = __a_get_dim(a).square();
            T* ad = updated(a);
            for(size_t i = 0; i < size; ++i) ad[i] = value; // not a memset due to complex
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct round_square : public kernel< round_square<T> > 
    {
        typedef void (round_square::*F)(const matrix<T>&, std::vector<T>*&);

        inline void l(const matrix<T>& a, std::vector<T>*& ac){
            pin(current(a));
        }

        inline void c(const matrix<T>& a, std::vector<T>*& ac){
            // gs
            __A_TIME_C("ambient_round_square_c_kernel"); 
            T* ad = current(a);
            size_t sizey = __a_get_dim(a).y;
            for(int i=0; i < sizey; i++){
                double v = std::abs(ad[i]);
                if(v > 1e-10) ac->push_back(v*v);
            }
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct cast_to_vector : public kernel< cast_to_vector<T> > 
    {
        typedef void (cast_to_vector::*F)(std::vector<T>*&, const matrix<T>&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(std::vector<T>*& ac, const matrix<T>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
            pin(current(a));
        }

        inline void c(std::vector<T>*& ac, const matrix<T>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
            // gs
            __A_TIME_C("ambient_cast_to_vector_c_kernel");
            T* ad = c_current(a);
            for(int j=0; j < n; ++j) memcpy((void*)&(*ac)[j*lda + offset],(void*)&ad[j*m], m*sizeof(T));  
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct cast_from_vector : public kernel< cast_from_vector<T> > 
    {
        typedef void (cast_from_vector::*F)(const std::vector<T>*&, matrix<T>&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(const std::vector<T>*& ac, matrix<T>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
            pin(current(a));
        }

        inline void c(const std::vector<T>*& ac, matrix<T>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
            __A_TIME_C("ambient_cast_from_vector_c_kernel"); 
            T* ad = updated(a);
            for(int j=0; j < n; ++j) memcpy((void*)&ad[j*m],(void*)&(*ac)[offset + j*lda], m*sizeof(T));
            __A_TIME_C_STOP 
        }
    };

    template<typename T1, typename T2>
    struct cast_from_vector_t : public kernel< cast_from_vector_t<T1,T2> > 
    {
        typedef void (cast_from_vector_t::*F)(const std::vector<T1>*&, matrix<T2>&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(const std::vector<T1>*& ac, matrix<T2>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
            pin(current(a));
        }

        inline void c(const std::vector<T1>*& ac, matrix<T2>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
            __A_TIME_C("ambient_cast_from_vector_t_c_kernel"); 
            T2* ad = updated(a);
            const T1* sd = &(*ac)[offset];
            for(int j=0; j < n; ++j) 
                for(int i=0; i < m; ++i)
                    ad[j*m + i] = sd[j*lda + i];
            __A_TIME_C_STOP 
        }
    };

    template<typename T>
    struct save : public kernel< save<T> >
    {
        typedef void (save::*F)(matrix<T> const&, const size_t& );

        inline void l(matrix<T> const& a, const size_t& tag){
            pin(current(a));
        }

        inline void c(matrix<T> const& a, const size_t& tag ){
            T* ad = (T*)current(a);
            ambient::io_manager.save(ad, tag, __a_sizeof(a)); 
        }
    };

    template<typename T>
    struct load : public kernel< load<T> >
    {
        typedef void (load::*F)(matrix<T>&, const size_t&);

        inline void l(matrix<T>& a, const size_t& tag){
            pin(current(a));
        }

        inline void c(matrix<T>& a, const size_t& tag){
            T* ad = (T*)updated(a);
            ambient::io_manager.load(ad, tag, __a_sizeof(a)); 
        }
    };

    template<typename T>
    struct validation : public kernel< validation<T> > 
    {
        typedef void (validation::*F)(const matrix<T>&, const matrix<T>&, future<bool>&);

        inline double distance(const std::complex<double>& a, const std::complex<double>& b){
            return fabs(std::norm(a) - std::norm(b));
        }
        inline double magnitude(const std::complex<double>& a, const std::complex<double>& b){
            return std::max(fabs(std::norm(a)), fabs(std::norm(b)));
        }
        inline double distance(double a, double b) { return fabs(fabs(a) - fabs(b));    }
        inline double magnitude(double a, double b){ return std::max(fabs(a), fabs(b)); }

        inline void l(const matrix<T>& a, const matrix<T>& b, future<bool>& ret){
            pin(current(a));  //if(!ctxt.involved()) return;
            assign(current(b)); 
        }

        inline void c(const matrix<T>& a, const matrix<T>& b, future<bool>& ret){ // see paper for Reference Dongara 
            __A_TIME_C("ambient_validation_kernel"); 
            T* ad = c_current(a); 
            T* bd = c_current(b); 
            double epsilon = std::numeric_limits<double>::epsilon();
            int count = 0;
            size_t sizey = std::min(__a_get_dim(a).y, __a_get_dim(b).y);
            size_t sizex = std::min(__a_get_dim(a).x, __a_get_dim(b).x);

            for(size_t i=0; i < sizey; ++i){
                for(size_t j=0; j < sizex; ++j){
                    T av = ad[i+j*__a_get_dim(a).y];
                    T bv = bd[i+j*__a_get_dim(b).y];
                    double d = distance(av, bv);
                    double m = magnitude(av, bv);
                    if(d > epsilon*256 && d/m > epsilon*256){ // || av*bv < 0 // 16 is recommended, 256 because MKL isn't bitwise stable
                        std::cout << i << " " << j << " : " << av << " " << bv << ", eps: " << d << std::endl;
                        ret.get_value() = false;
                        if(++count > 10) return;
                    }
                }
            }
            __A_TIME_C_STOP 
        }
    };

    template<typename T>
    struct labrd_update_col : public kernel< labrd_update_col<T> > 
    { // svd
        typedef void (labrd_update_col::*F)(matrix<T>& , const matrix<T>&, matrix<T>&, const matrix<T>&, matrix<T>&, matrix<T>&, const int&);
        inline void l(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, matrix<T>& tq, matrix<T>& d, const int& i){
            pin(current(say));
            assign(current(sax));
            assign(current(sy));
            assign(current(sx));
            assign(current(tq));
            assign(current(d));
        }
        inline void c(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, matrix<T>& tq, matrix<T>& d, const int& i){
            __A_TIME_C("ambient_labrd_update_col_c_kernel"); 
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

                __a_memptf_r<T, __a_memcpy>(sayd, ldsay, dim2(i, i-1), 
                                            saxd, ldsax, dim2(i, i-1), 
                                            dim2( num_cols(say)-i, 1));

                dgemv_("N", &ri, &i, &mone, &sayd[ i ], &ldsay, &syd[ i ], &ldsy, &one, &sayd[i + i*ldsay], &lone);
                dgemv_("N", &ri, &i, &mone, &sxd[ i ], &ldsx, &sayd[ i*ldsay ], &lone, &one, &sayd[i + i*ldsay], &lone);
                
                dlarfg_( &ri, &sayd[i+i*ldsay], &sayd[std::min(i+1, m-1)+i*ldsay], &lone, &tqd[i] );
                dd[i] = sayd[i+i*ldsay];
                sayd[i+i*ldsay] = 1.000;
                
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct labrd_reduce_col : public kernel< labrd_reduce_col<T> > 
    { // svd
        typedef void (labrd_reduce_col::*F)(matrix<T>& , const matrix<T>&, matrix<T>&, const matrix<T>&, const int&);
        inline void l(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, const int& i){
            pin(current(say));
            assign(current(sax));
            assign(current(sy));
            assign(current(sx));
        }
        inline void c(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, const int& i){
            __A_TIME_C("ambient_labrd_reduce_col_c_kernel"); 
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

            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct labrd_update_row : public kernel< labrd_update_row<T> > 
    { // svd
        typedef void (labrd_update_row::*F)(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, matrix<T>& tp, matrix<T>& e, const int& i);
        inline void l(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, matrix<T>& tp, matrix<T>& e, const int& i){
            pin(current(say));
            assign(current(sax));
            assign(current(sy));
            assign(current(sx));
            assign(current(tp));
            assign(current(e));
        }
        inline void c(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, matrix<T>& tp, matrix<T>& e, const int& i){
            __A_TIME_C("ambient_labrd_update_row_c_kernel"); 
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
                
                __a_memptf_r<T, __a_memcpy>(saxd, ldsax, dim2(i, i), 
                                            sayd, ldsay, dim2(i, i), 
                                            dim2( 1, ldsax-i ));
                
                dgemv_("T", &i, &rj, &mone, &saxd[(i+1)*ldsax], &ldsax, &sxd[i], &ldsx, &one, &saxd[ i + (i+1)*ldsax], &ldsax);
                dgemv_("N", &rj, &r3, &mone, &syd[ i+1 ], &ldsy, &saxd[i], &ldsax, &one, &saxd[i + (i+1)*ldsax], &ldsax);

                dlarfg_(&rj, &saxd[i + (i+1)*ldsax], &saxd[i + std::min(i+2, n-1)*ldsax], &ldsax, &tpd[i] );
                ed[i] = saxd[i + (i+1)*ldsax];
                saxd[i + (i+1)*ldsax] = 1.000;

            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct labrd_reduce_row : public kernel< labrd_reduce_row<T> > 
    { // svd
        typedef void (labrd_reduce_row::*F)(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, const int& i);
        inline void l(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, const int& i){
            pin(current(say));
            assign(current(sax));
            assign(current(sy));
            assign(current(sx));
        }
        inline void c(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, const int& i){
            __A_TIME_C("ambient_labrd_reduce_row_c_kernel"); 
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
            __A_TIME_C_STOP
        }
    };

    // {{{ MKL LAPACK kernels

    template<typename T, PLASMA_enum TR>
    struct larfg : public kernel< larfg<T,TR> > 
    { // svd
        typedef void (larfg::*F)(matrix<T>&, matrix<T>&, matrix<T>&, const size_t&);
        inline void l(matrix<T>& a, matrix<T>& t, matrix<T>& d, const size_t& k){
            pin(current(a));
            assign(current(t));
            assign(current(d));
        }
        inline void c(matrix<T>& a, matrix<T>& t, matrix<T>& d, const size_t& k){
            __A_TIME_C("ambient_larfg_c_kernel"); 
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
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct gebd2 : public kernel< gebd2<T> > 
    { // gs
        typedef void (gebd2::*F)(matrix<T>&, matrix<T>&, matrix<T>&, matrix<T>&, matrix<T>&);
        inline void l(matrix<T>& a, matrix<T>& d, matrix<T>& e, matrix<T>& tq, matrix<T>& tp){
            pin(current(a));
            assign(current(d));
            assign(current(e));
            assign(current(tq));
            assign(current(tp));
        }
        inline void c(matrix<T>& a, matrix<T>& d, matrix<T>& e, matrix<T>& tq, matrix<T>& tp){
            __A_TIME_C("ambient_gebd2_c_kernel"); 
            int m = a.num_rows();
            int n = a.num_cols();
            int lda = a.num_rows();
            int info;

            T* work = (T*)malloc(std::max(m,n)*sizeof(T));
            dgebd2_(&m, &n, (T*)s_updated(a), &lda, (T*)updated(d), (T*)updated(e), (T*)updated(tq), (T*)updated(tp), work, &info);
            free(work);
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct gebrd : public kernel< gebrd<T> > 
    { // gs
        typedef void (gebrd::*F)(matrix<T>&, weak_view<T>&, weak_view<T>&, weak_view<T>&, weak_view<T>&);
        inline void l(matrix<T>& a, weak_view<T>& d, weak_view<T>& e, weak_view<T>& q, weak_view<T>& p){
            pin(current(a));
            assign(current(d));
            assign(current(e));
            assign(current(q));
            assign(current(p));
        }
        inline void c(matrix<T>& a, weak_view<T>& d, weak_view<T>& e, weak_view<T>& q, weak_view<T>& p){
            __A_TIME_C("ambient_gebrd_c_kernel");
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
            memcpy(ac, ad, m*n*sizeof(T));

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
                memcpy(&pd[n*j], &ac[lda*j], sizeof(T)*k);
            }

            free(work); free(ac);
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct gbbrd : public kernel< gbbrd<T> > 
    { // gs
        typedef void (gbbrd::*F)(matrix<T>&, weak_view<T>&, weak_view<T>&, weak_view<T>&, weak_view<T>&);
        inline void l(matrix<T>& a, weak_view<T>& d, weak_view<T>& e, weak_view<T>& q, weak_view<T>& p){
            pin(current(a));
            assign(current(d));
            assign(current(e));
            assign(current(q));
            assign(current(p));
        }
        inline void c(matrix<T>& a, weak_view<T>& d, weak_view<T>& e, weak_view<T>& q, weak_view<T>& p){
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

            __A_TIME_C("ambient_gbbrd_c_kernel"); 
            dgbbrd_("B", &m, &n, &zero, &kl, &ku, ad, &k, dd, ed, 
                    qd, &m, pd, &n, NULL, &one, work, &info);

            __A_TIME_C_STOP
            free(work);
        }
    };

    template<typename T>
    struct bdsqr : public kernel< bdsqr<T> > 
    { // gs
        typedef void (bdsqr::*F)(matrix<T>&, matrix<T>&, matrix<T>&, matrix<T>&);
        inline void l(matrix<T>& d, matrix<T>& e, matrix<T>& u, matrix<T>& v){
            pin(current(d));
            assign(current(e));
            assign(current(u));
            assign(current(v));
        }
        inline void c(matrix<T>& d, matrix<T>& e, matrix<T>& u, matrix<T>& v){
            __A_TIME_C("ambient_bdsqr_c_kernel"); 
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

            __A_TIME_C_STOP
        }
    };
    
    template<typename T>
    struct svd : public kernel< svd<T> > 
    {
        typedef void (svd::*F)(const matrix<T>&, weak_view<T>&, weak_view<T>&, weak_view<double>&);

        inline void l(const matrix<T>& a, weak_view<T>& u, weak_view<T>& vt, weak_view<double>& s){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(s));
            assign(current(u));
            assign(current(vt));
        }

        inline void c(const matrix<T>& a, weak_view<T>& u, weak_view<T>& vt, weak_view<double>& s){
            // gs
            __A_TIME_C("ambient_svd_c_kernel"); 
            int m = __a_get_dim(a).y;
            int n = __a_get_dim(a).x;
            int k = std::min(m,n);
            int info;
            int lwork = -1; // C - Alex, netlib said -1 for the best workspace
            T wkopt;
            T* ad  = current(a);
            T* ud  = updated(u);
            T* vtd = updated(vt);
            double* sd  = updated(s);
            double* rwork; // = new double[5*m]; // C - useless for double but need for complex 
            T* work;
            gesvd( "S", "S", &m, &n, ad, &m, sd, ud, &m, vtd, &k, &wkopt, &lwork, rwork, &info ); // query and allocate the optimal workspace
            lwork = OptimalSize(wkopt);
            work = (T*)malloc( lwork*sizeof(T) );
            gesvd( "S", "S", &m, &n, ad, &m, sd, ud, &m, vtd, &k, work, &lwork, rwork, &info );   // compute SVD
            assert( info == 0 ); // otherwise the algorithm computing atomic SVD failed to converge
            free(work);
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct qr : public kernel< qr<T> > 
    {
        typedef void (qr::*F)(const matrix<T>&, weak_view<T>&, weak_view<T>&);

        inline void l(const matrix<T>& a, weak_view<T>& q, weak_view<T>& r){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(q));
            assign(current(r));
        }

        inline void c(const matrix<T>& a, weak_view<T>& q, weak_view<T>& r){
            // gs
            __A_TIME_C("ambient_qr_c_kernel"); 
            int m = __a_get_dim(a).y; //numrow a
            int n = __a_get_dim(a).x; //numcol a, numcol r
            int k = std::min(m,n); //numrow r
            int info;
            int lwork = -1; 
            T wkopt;
            T* tau = (T*)malloc(k*sizeof(T));
            T* ad  = current(a);
            T* qd  = updated(q);
            T* rd = p_updated(r);
            T* work;
            T* more_work;
            T  kwork;

            geqrf(&m, &n, ad, &m, tau, &kwork, &lwork, &info);
            lwork = OptimalSize(kwork);
            work = (T*)malloc( lwork*sizeof(T) );
            geqrf(&m, &n, ad, &m, tau, work, &lwork, &info);
            assert( info == 0 );

            for (std::size_t c = 0; c < n; ++c)
                for (std::size_t r = 0; r <= c && r < k; ++r)
                    rd[r+k*c] = ad[r+m*c]; 

            lwork = -1;

            getq_qr(&m, &k, &k, ad, &m, tau, &kwork, &lwork, &info);

            lwork = OptimalSize(kwork);
            more_work = (T*)malloc( lwork*sizeof(T) );
            getq_qr(&m, &k, &k, ad, &m, tau, more_work, &lwork, &info);
            assert( info == 0 ); 
             
            memcpy((void*)qd, (void*)ad, k*__a_get_dim(a).y*sizeof(T)); // l 235 

            free(work);
            free(more_work);
            free(tau);

            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct lq : public kernel< lq<T> > 
    {
        typedef void (lq::*F)(const matrix<T>&, weak_view<T>&, weak_view<T>&);

        inline void l(const matrix<T>& a, weak_view<T>& l, weak_view<T>& q){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(l));
            assign(current(q));
        }

        inline void c(const matrix<T>& a, weak_view<T>& l, weak_view<T>& q){
            // gs
            __A_TIME_C("ambient_lq_c_kernel"); 
            int m = __a_get_dim(a).y; //numrow a, numrow l
            int n = __a_get_dim(a).x; //numcol a
            int k = std::min(m,n); //numcol l
            int info;
            int lwork = -1; 
            T wkopt;
            T* tau = (T*)malloc(k*sizeof(T));
            T* ad  = current(a);
            T* ld  = p_updated(l);
            T* qd  = updated(q);
            T* work;
            T* more_work;
            T  kwork;

            gelqf(&m, &n, ad, &m, tau, &kwork, &lwork, &info);
            lwork = OptimalSize(kwork);
            work = (T*)malloc( lwork*sizeof(T) );
            gelqf(&m, &n, ad, &m, tau, work, &lwork, &info);
            assert( info == 0 );

            for (std::size_t c = 0; c < k; ++c)
                for (std::size_t r = c; r < m ;++r)
                    ld[r+m*c] = ad[r+m*c]; 

            lwork = -1;
            getq_lq(&k, &n, &k, ad, &m, tau, &kwork, &lwork, &info);

            lwork = OptimalSize(kwork);
            more_work = (T*)malloc( lwork*sizeof(T) );

            getq_lq(&k, &n, &k, ad, &m, tau, more_work, &lwork, &info);
            assert( info == 0 ); 

            for (std::size_t c = 0; c < n; ++c)
                for (std::size_t r = 0; r < k ;++r)
                    qd[r+k*c] = ad[r+m*c]; 

            free(work);
            free(more_work);
            free(tau);

            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct heev : public kernel< heev<T> > 
    {
        typedef void (heev::*F)(matrix<T>&, weak_view<double>&);

        inline void l(matrix<T>& a, weak_view<double>& w){
            pin(current(a)); //if(!ctxt.involved()) return;
            assign(current(w));
        }

        inline void c(matrix<T>& a, weak_view<double>& w){
            // gs
            __A_TIME_C("ambient_heev_c_kernel"); 
            int m = __a_get_dim(a).y;
            int info, lwork = -1;
            double wkopt;
            double* work;
            double* ad = (double*)__a_solidify(c_current(a), __a_sizeof(a));
            double* wd = (double*)__a_solidify(c_current(w), __a_sizeof(w));

            dsyev_("V","U",&m,ad,&m,wd,&wkopt,&lwork,&info);
            lwork = (int)wkopt;
            work = (double*)malloc( lwork*sizeof(double) );
            dsyev_("V","U",&m,ad,&m,wd,work,&lwork,&info);
            assert( info == 0 ); // otherwise the algorithm computing SYEV failed to converge
            // First we reverse the eigenvalues, to be coherent with the serial version ! 
            for(int i=0; i < (int)(m/2); i++){ 
                wkopt = wd[i];
                wd[i] = wd[m-i-1];
                wd[m-i-1] = wkopt;
            } 
            // Second we reverse the eigenvectors
            size_t len = m*sizeof(double);
            work = (double*)realloc(work, len);
            for (int i=0; i < (int)(m/2); i++){ 
                memcpy(work, &ad[i*m], len);
                memcpy(&ad[i*m], &ad[(m-1-i)*m], len);
                memcpy(&ad[(m-1-i)*m], work, len);
            }
            __a_disperse(ad, updated(a), __a_sizeof(a));
            __a_disperse(wd, updated(w), __a_sizeof(w));
            free(work);
            __A_TIME_C_STOP
        }
    };

    // }}}
} } }

#undef pin
#endif
