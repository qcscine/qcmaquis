#ifndef __AMBIENT_NUMERIC_MATRIX_KERNELS_ATOMICS_HPP__
#define __AMBIENT_NUMERIC_MATRIX_KERNELS_ATOMICS_HPP__

extern "C" {
    void dgemm_(const char*,const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
}

namespace ambient { namespace numeric { namespace kernels {

    using ambient::numeric::matrix_impl;

    template<class Tag1, class Tag2, typename T>
    struct gemm_general_atomic : public kernel_atomic< gemm_general_atomic<Tag1, Tag2, T> > 
    { // gs
        typedef void(gemm_general_atomic::*F)(const matrix_impl<T>&, const matrix_impl<T>&, matrix_impl<T>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<T>& b, matrix_impl<T>& c){
            this->ctxt_select("1 from ambient as gemm_general_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
            this->assign(ui_l_current(c));
        }
        inline void c(const matrix_impl<double>& a, const matrix_impl<double>& b, matrix_impl<double>& c){
            __A_TIME_C("ambient_gemm_general_atomic_c_kernel"); 
            double* ad = ui_c_current(a)(0,0);
            double* bd = ui_c_current(b)(0,0);
            double* cd = ui_w_updated(c)(0,0);
            int m = Tag1::first(a);
            int n = Tag2::second(b);
            int k = Tag1::second(a);
            int lda = ui_c_get_dim(a).y;
            int ldb = ui_c_get_dim(b).y;
            int ldc = ui_c_get_dim(c).y;
            static const double alpha(1.0); 
            static const double beta(0.0);
            dgemm_(Tag1::code(), Tag2::code(), &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
            __A_TIME_C_STOP
        }
        inline void c(const matrix_impl<std::complex<double> >& a, const matrix_impl<std::complex<double> >& b, matrix_impl<std::complex<double> >& c){
            __A_TIME_C("ambient_gemm_general_atomic_c_kernel"); 
            T* ad   = ui_c_current(a)(0,0);
            T* bd   = ui_c_current(b)(0,0);
            T* cd   = ui_w_updated(c)(0,0);
            int m   = ui_c_get_dim(a).y;
            int n   = ui_c_get_dim(b).x;
            int k   = ui_c_get_dim(b).y;
            T alpha(1.0); 
            T beta(0.0);
            gemm("N","N", &m, &n, &k, &alpha, ad, &m, bd, &k, &beta, cd, &m);
            __A_TIME_C_STOP
        }
    };


    template<typename T>
    struct copy_atomic : public kernel_atomic< copy_atomic<T> > 
    { // gs
        typedef void(copy_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&);

        inline void l(matrix_impl<T>& ac, const matrix_impl<T>& a){
            this->ctxt_select("1 from ambient as copy_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(ac));
        }

        inline void c(matrix_impl<T>& ac, const matrix_impl<T>& a){
            __A_TIME_C("ambient_copy_atomic_c_kernel"); 
            T* ad  = ui_c_current(a)(0,0);
            T* acd  = ui_w_updated(ac)(0,0);
            memcpy(acd, ad, ui_c_get_mem_size(a));
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct reshape_l2r_atomic : public kernel_atomic< reshape_l2r_atomic<T> > 
    { // gs
        typedef void (reshape_l2r_atomic::*F)(const matrix_impl<T>&, matrix_impl<T>&,
                                              const size_t&, const size_t&, 
                                              const size_t&, const size_t&, const size_t&);

        inline void l(const matrix_impl<T>& left, matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as reshape_l2r_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(right));
            this->assign(ui_l_current(left));
        }

        inline void c(const matrix_impl<T>& left, matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            __A_TIME_C("ambient_reshape_l2r_atomic_c_kernel"); 
            __a_atomic_refresh(right); // refreshing updated memory
            for(size_t ss = 0; ss < sdim; ++ss){
                __a_memptf_atomic_r<T, __a_memcpy>(right, dim2(ss*rdim + right_offset, 0), 
                                                   left,  dim2(0, ss*ldim + left_offset), 
                                                   dim2( rdim, ldim ));
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct reshape_r2l_atomic : public kernel_atomic< reshape_r2l_atomic<T> > 
    { // gs
        typedef void (reshape_r2l_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&,
                                              const size_t&, const size_t&, 
                                              const size_t&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& left, const matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as reshape_l2r_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(left));
            this->assign(ui_l_current(right));
        }

        inline void c(matrix_impl<T>& left, const matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            __A_TIME_C("ambient_reshape_r2l_atomic_c_kernel"); 
            __a_atomic_refresh(left); // refreshing updated memory
            for(size_t ss = 0; ss < sdim; ++ss)
                __a_memptf_atomic_r<T, __a_memcpy>(left,  dim2(0, ss*ldim + left_offset), 
                                                   right, dim2(ss*rdim + right_offset,0), 
                                                   dim2( rdim, ldim ));
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct lb_tensor_mpo_atomic : public kernel_atomic< lb_tensor_mpo_atomic<T> > 
    { // gs
        typedef void (lb_tensor_mpo_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&, const matrix_impl<T>&,
                                                const size_t&, const size_t&, 
                                                const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as rb_tensor_mpo_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(out));
            this->assign(ui_l_current(in));
            this->assign(ui_l_current(alfa));
        }

        inline void c(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            __A_TIME_C("ambient_lb_tensor_mpo_atomic_c_kernel"); 
            __a_atomic_refresh(out); // refreshing updated memory
            T* alfad = ui_c_current(alfa)(0,0);
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T alfa_t = alfad[ss1 + ui_c_get_dim(alfa).y*ss2];
                    __a_memptf_atomic_r<T, __a_memscal>(out, dim2(0, out_offset + ss2*ldim),
                                                        in,  dim2(0, in_offset + ss1*ldim),
                                                        dim2(rdim, ldim), alfa_t);
                }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct rb_tensor_mpo_atomic : public kernel_atomic< rb_tensor_mpo_atomic<T> > 
    { // gs
        typedef void (rb_tensor_mpo_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&, const matrix_impl<T>&,
                                                const size_t&, const size_t&, 
                                                const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as rb_tensor_mpo_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(out));
            this->assign(ui_l_current(in));
            this->assign(ui_l_current(alfa));
        }

        inline void c(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            __A_TIME_C("ambient_rb_tensor_mpo_atomic_c_kernel"); 
            __a_atomic_refresh(out); // refreshing updated memory
            T* alfad = ui_c_current(alfa)(0,0);
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T alfa_t = alfad[ss1 + ui_c_get_dim(alfa).y*ss2];
                    __a_memptf_atomic_r<T, __a_memscal>(out, dim2(out_offset + ss2*rdim, 0),
                                                        in,  dim2(in_offset + ss1*rdim, 0),
                                                        dim2(rdim, ldim), alfa_t);
                }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scalar_norm_atomic : public kernel_atomic< scalar_norm_atomic<T> > 
    {// gs
        typedef void (scalar_norm_atomic::*F)(const matrix_impl<T>&, future<double>&);

        inline void l(const matrix_impl<T>& a, future<double>& norm){
            this->ctxt_select("* from ambient as scalar_norm_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(const matrix_impl<T>& a, future<double>& norm){
            __A_TIME_C("ambient_scalar_norm_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            norm.get_value() = alps::numeric::real(__a_dot(ad, ad, ui_c_get_dim(a).square()));
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct overlap_atomic : public kernel_atomic< overlap_atomic<T> > 
    { // gs
        typedef void (overlap_atomic::*F)(const matrix_impl<T>&, const matrix_impl<T>&, future<T>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<T>& b, future<T>& overlap){
            this->ctxt_select("* from ambient as overlap_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
        }

        inline void c(const matrix_impl<T>& a, const matrix_impl<T>& b, future<T>& overlap){
            __A_TIME_C("ambient_scalar_overlap_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* bd = ui_c_current(b)(0,0);
            overlap.get_value() = __a_dot(ad, bd, ui_c_get_dim(a).square());
            __A_TIME_C_STOP
        }
    };

        
    template<typename T>
    struct add_atomic : public kernel_atomic< add_atomic<T> > 
    { // gs
        typedef void (add_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&);

        inline void l(matrix_impl<T>& a, const matrix_impl<T>& b){
            this->ctxt_select("1 from ambient as add_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
        }

        inline void c(matrix_impl<T>& a, const matrix_impl<T>& b){
            __A_TIME_C("ambient_add_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* bd = ui_c_current(b)(0,0);
            T* ar = ui_r_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            static const T sign = 1.;
            static const int ONE = 1;
            axpy(&size, &sign, bd, &ONE, ar, &ONE);
            __A_TIME_C_STOP
        }
    };

        
    template<typename T>
    struct sub_atomic : public kernel_atomic< sub_atomic<T> > 
    { // gs
        typedef void (sub_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&);

        inline void l(matrix_impl<T>& a, const matrix_impl<T>& b){
            this->ctxt_select("1 from ambient as sub_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
        }

        inline void c(matrix_impl<T>& a, const matrix_impl<T>& b){
            __A_TIME_C("ambient_sub_atomic_c_kernel"); 
            T* bd = ui_c_current(b)(0,0);
            T* ar = ui_r_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            static const T sign = -1.;
            static const int ONE = 1;
            axpy(&size, &sign, bd, &ONE, ar, &ONE);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scale_atomic : public kernel_atomic< scale_atomic<T> > 
    { // gs
        typedef void (scale_atomic::*F)(matrix_impl<T>&, const future<T>&);

        inline void l(matrix_impl<T>& a, const future<T>& t){
            this->ctxt_select("1 from ambient as scale_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(matrix_impl<double>& a, const future<double>& t){
            __A_TIME_C("ambient_scale_atomic_c_kernel"); 
            T* ar = ui_r_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            static const int ONE = 1;
            dscal_( &size, &t.get_value(), ar, &ONE );
            __A_TIME_C_STOP
        }

        inline void c(matrix_impl<std::complex<double> >& a, const future< std::complex<double> >& t){
            __A_TIME_C("ambient_scale_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* ar = ui_w_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            for(int k=0; k < size; k++) 
                ar[k] = ad[k] * t.get_value();
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scale_inverse_atomic : public kernel_atomic< scale_inverse_atomic<T> > 
    { // gs
        typedef void (scale_inverse_atomic::*F)(matrix_impl<T>&, const future<T>&);

        inline void l(matrix_impl<T>& a, const future<T>& t){
            this->ctxt_select("1 from ambient as scale_inverse_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(matrix_impl<double>& a, const future<double>& t){
            __A_TIME_C("ambient_scale_inverse_atomic_c_kernel"); 
            T* ar = ui_r_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            static const int ONE = 1;
            double factor = 1. / t.get_value();
            dscal_( &size, &factor, ar, &ONE );
            __A_TIME_C_STOP
        }

        inline void c(matrix_impl<std::complex<double> >& a, const future< std::complex<double> >& t){
            __A_TIME_C("ambient_scale_inverse_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* ar = ui_w_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            for(int k=0; k < size; k++) 
                ar[k] = ad[k] / t.get_value();
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct transpose_out_atomic : public kernel_atomic< transpose_out_atomic<T> > 
    { // gs
        typedef void (transpose_out_atomic::*F)(const matrix_impl<T>&, matrix_impl<T>&);

        inline void l(const matrix_impl<T>& a, matrix_impl<T>& t){
            this->ctxt_select("1 from ambient as transpose_out_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(t));
        }

        inline void c(const matrix_impl<T>& a, matrix_impl<T>& t){
            __A_TIME_C("ambient_transpose_out_atomic_c_kernel"); 
            T* od = ui_c_current(a)(0,0);
            T* td = ui_w_updated(t)(0,0);
            int m = ui_c_get_dim(a).y;
            int n = ui_c_get_dim(a).x;

            for(int i = 0; i < m; i++){
                for(int j = 0; j < n; j++) *td++ = od[j*m];
                od++;
            }
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct resize_atomic : public kernel_atomic< resize_atomic<T> > 
    {
        typedef void (resize_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& r, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            this->ctxt_select("1 from ambient as resize_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(r));
        }

        inline void c(matrix_impl<T>& r, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            __A_TIME_C("ambient_resize_atomic_c_kernel"); 
            __a_memptf_atomic_r<T, __a_memcpy>(r, dim2(0,0), a, dim2(0,0), dim2(n, m)); 
            __A_TIME_C_STOP
        }
    };
        
    // {{{ MKL LAPACK kernels

    template<typename T>
    struct svd_atomic : public kernel_unpinned< svd_atomic<T> > 
    {
        typedef void (svd_atomic::*F)(const matrix_impl<T>&, matrix_impl<T>&, matrix_impl<T>&, matrix_impl<double>&);

        inline void l(const matrix_impl<T>& a, matrix_impl<T>& u, matrix_impl<T>& vt, matrix_impl<double>& s)
        {
            this->ctxt_select("1 from ambient as svd_atomic"); //if(!ctxt.involved()) return;
            this->atomic_conditional_assign(ui_l_current(s));
            this->atomic_conditional_assign(ui_l_current(a));
            this->atomic_conditional_assign(ui_l_current(u));
            this->atomic_conditional_assign(ui_l_current(vt));
        }

        inline void c(const matrix_impl<T>& a, matrix_impl<T>& u, matrix_impl<T>& vt, matrix_impl<double>& s)
        { // gs
            __A_TIME_C("ambient_svd_atomic_c_kernel"); 
            int m = ui_c_get_dim(a).y;
            int n = ui_c_get_dim(a).x;
            int k = std::min(m,n);
            int info;
            int lwork = -1; // C - Alex, netlib said -1 for the best workspace
            T wkopt;
            T* ad  = ui_c_current(a)(0,0);
            T* ud  = ui_r_updated(u)(0,0);
            T* vtd = ui_r_updated(vt)(0,0);
            double* sd  = ui_r_updated(s)(0,0);
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
    struct heev_atomic : public kernel_unpinned< heev_atomic<T> > 
    {
        typedef void (heev_atomic::*F)(matrix_impl<T>&, matrix_impl<double>&);

        inline void l(matrix_impl<T>& a, matrix_impl<double>& w){
            this->ctxt_select("1 from ambient as heev_atomic"); //if(!ctxt.involved()) return;
            this->atomic_conditional_assign(ui_l_current(a));
            this->atomic_conditional_assign(ui_l_current(w));
        }

        inline void c(matrix_impl<T>& a, matrix_impl<double>& w){
            // gs
            __A_TIME_C("ambient_heev_atomic_c_kernel"); 
            int m = ui_c_get_dim(a).y;
            int info, lwork = -1;
            double wkopt;
            double* work;
            double* ad = (double*)__a_solidify_atomic<T>(a);
            double* wd = (double*)__a_solidify_atomic<T>(w);

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
            __a_disperse_atomic<T>(ad, a);
            __a_disperse_atomic<T>(wd, w);
            free(work);
            __A_TIME_C_STOP
        }
    };

    // }}}
} } }
#endif
