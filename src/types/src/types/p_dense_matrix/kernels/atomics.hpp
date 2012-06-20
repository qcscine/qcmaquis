#ifndef __MAQUIS_TYPES_KERNELS_ATOMICS_HPP__
#define __MAQUIS_TYPES_KERNELS_ATOMICS_HPP__

extern "C" {
    void dgemm_(const char*,const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
}


namespace ambient {

    template<typename T>
    struct gemm_general_atomic : public ambient::kernel_atomic< gemm_general_atomic<T> > 
    { // gs
        typedef void(gemm_general_atomic::*F)(const maquis::types::p_dense_matrix_impl<T>&, 
                         const maquis::types::p_dense_matrix_impl<T>&, 
                               maquis::types::p_dense_matrix_impl<T>&);

        inline void l(const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
            this->ctxt_select("1 from ambient as gemm_general_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
            this->assign(ui_l_current(c));
        }
        inline void c(const maquis::types::p_dense_matrix_impl<double>& a, const maquis::types::p_dense_matrix_impl<double>& b, maquis::types::p_dense_matrix_impl<double>& c){
            __A_TIME_C("ambient_gemm_general_atomic_c_kernel"); 
            double* ad = ui_c_current(a)(0,0);
            double* bd = ui_c_current(b)(0,0);
            double* cd = ui_w_updated(c)(0,0);
            int m = ui_c_get_dim(a).y;
            int n = ui_c_get_dim(b).x;
            int k = ui_c_get_dim(b).y;
            static const double alpha(1.0); 
            static const double beta(0.0);
            dgemm_("N","N", &m, &n, &k, &alpha, ad, &m, bd, &k, &beta, cd, &m);
            __A_TIME_C_STOP
        }
        inline void c(const maquis::types::p_dense_matrix_impl<std::complex<double> >& a, const maquis::types::p_dense_matrix_impl<std::complex<double> >& b, maquis::types::p_dense_matrix_impl<std::complex<double> >& c){
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
    struct copy_atomic : public ambient::kernel_atomic< copy_atomic<T> > 
    { // gs
        typedef void(copy_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& ac, const maquis::types::p_dense_matrix_impl<T>& a){
            this->ctxt_select("1 from ambient as copy_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(ac));
        }

        inline void c(maquis::types::p_dense_matrix_impl<T>& ac, const maquis::types::p_dense_matrix_impl<T>& a){
            __A_TIME_C("ambient_copy_atomic_c_kernel"); 
            T* ad  = ui_c_current(a)(0,0);
            T* acd  = ui_w_updated(ac)(0,0);
            memcpy(acd, ad, ui_c_get_mem_size(a));
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct reshape_l2r_atomic : public ambient::kernel_atomic< reshape_l2r_atomic<T> > 
    { // gs
        typedef void (reshape_l2r_atomic::*F)(const maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(const maquis::types::p_dense_matrix_impl<T>& left, maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as reshape_l2r_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(right));
            this->assign(ui_l_current(left));
        }

        inline void c(const maquis::types::p_dense_matrix_impl<T>& left, maquis::types::p_dense_matrix_impl<T>& right,
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
    struct reshape_r2l_atomic : public ambient::kernel_atomic< reshape_r2l_atomic<T> > 
    { // gs
        typedef void (reshape_r2l_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
               const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& left, const maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as reshape_l2r_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(left));
            this->assign(ui_l_current(right));
        }

        inline void c(maquis::types::p_dense_matrix_impl<T>& left, const maquis::types::p_dense_matrix_impl<T>& right,
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
    struct lb_tensor_mpo_atomic : public ambient::kernel_atomic< lb_tensor_mpo_atomic<T> > 
    { // gs
        typedef void (lb_tensor_mpo_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
               const size_t&, const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as rb_tensor_mpo_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(out));
            this->assign(ui_l_current(in));
            this->assign(ui_l_current(alfa));
        }

        inline void c(maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
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
    struct rb_tensor_mpo_atomic : public ambient::kernel_atomic< rb_tensor_mpo_atomic<T> > 
    { // gs
        typedef void (rb_tensor_mpo_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as rb_tensor_mpo_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(out));
            this->assign(ui_l_current(in));
            this->assign(ui_l_current(alfa));
        }

        inline void c(maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
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
    struct scalar_norm_atomic : public ambient::kernel_atomic< scalar_norm_atomic<T> > 
    {// gs
        typedef void (scalar_norm_atomic::*F)(const maquis::types::p_dense_matrix_impl<T>&, T*&);

        inline void l(const maquis::types::p_dense_matrix_impl<T>& a, T*& norm){
            this->ctxt_select("* from ambient as scalar_norm_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(const maquis::types::p_dense_matrix_impl<T>& a, T*& norm){
            __A_TIME_C("ambient_scalar_norm_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            *norm += __a_dot(ad, ad, ui_c_get_dim(a).square());
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scalar_overlap_atomic : public ambient::kernel_atomic< scalar_overlap_atomic<T> > 
    { // gs
        typedef void (scalar_overlap_atomic::*F)(const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, T*&);

        inline void l(const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, T*& overlap){
            this->ctxt_select("* from ambient as scalar_overlap_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
        }

        inline void c(const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, T*& overlap){
            __A_TIME_C("ambient_scalar_overlap_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* bd = ui_c_current(b)(0,0);
            *overlap += __a_dot(ad, bd, ui_c_get_dim(a).square());
            __A_TIME_C_STOP
        }
    };

        
    template<typename T>
    struct add_atomic : public ambient::kernel_atomic< add_atomic<T> > 
    { // gs
        typedef void (add_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            this->ctxt_select("1 from ambient as add_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
        }

        inline void c(maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            __A_TIME_C("ambient_add_atomic_c_kernel"); 
            T* bd = ui_c_current(b)(0,0);
            T* ar = ui_r_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            //if(ad != ar){
            //    for(int k = 0; k < size; k++)
            //        ar[k] = ad[k] + bd[k];
            //}else{
                static const T sign = 1.;
                static const int ONE = 1;
                axpy(&size, &sign, bd, &ONE, ar, &ONE);
            //}
            __A_TIME_C_STOP
        }
    };

        
    template<typename T>
    struct sub_atomic : public ambient::kernel_atomic< sub_atomic<T> > 
    { // gs
        typedef void (sub_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            this->ctxt_select("1 from ambient as sub_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
        }

        inline void c(maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            __A_TIME_C("ambient_sub_atomic_c_kernel"); 
            T* bd = ui_c_current(b)(0,0);
            T* ar = ui_r_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            //if(ad != ar){
            //    for(int k = 0; k < size; k++)
            //        ar[k] = ad[k] - bd[k];
            //}else{
                static const T sign = -1.;
                static const int ONE = 1;
                axpy(&size, &sign, bd, &ONE, ar, &ONE);
            //}
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scale_atomic : public ambient::kernel_atomic< scale_atomic<T> > 
    { // gs
        typedef void (scale_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, const T*&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& a, const T*& t){
            this->ctxt_select("1 from ambient as scale_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(maquis::types::p_dense_matrix_impl<double>& a, const double*& t){
            __A_TIME_C("ambient_scale_atomic_c_kernel"); 
            T* ar = ui_r_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            //if(ad != ar){
            //    for(int k=0; k < size; k++) 
            //        ar[k] = ad[k] * (*t);
            //}else{
                static const int ONE = 1;
                dscal_( &size, t, ar, &ONE );
            //}
            __A_TIME_C_STOP
        }

        inline void c(maquis::types::p_dense_matrix_impl<std::complex<double> >& a, const std::complex<double>*& t){
            __A_TIME_C("ambient_scale_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* ar = ui_w_updated(a)(0,0);
            int size = ui_c_get_dim(a).square();
            //if(ad != ar){
                for(int k=0; k < size; k++) 
                    ar[k] = ad[k] * (*t);
            //}else{
            //    int ONE = 1;
            //    zscal_( &size, t, ar, &ONE );
            //}
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct transpose_out_atomic : public ambient::kernel_atomic< transpose_out_atomic<T> > 
    { // gs
        typedef void (transpose_out_atomic::*F)(const maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<T>&);

        inline void l(const maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<T>& t){
            this->ctxt_select("1 from ambient as transpose_out_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(t));
        }

        inline void c(const maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<T>& t){
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
    struct resize_atomic : public ambient::kernel_atomic< resize_atomic<T> > 
    {
        typedef void (resize_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& r, const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
            this->ctxt_select("1 from ambient as resize_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(r));
        }

        inline void c(maquis::types::p_dense_matrix_impl<T>& r, const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
            __A_TIME_C("ambient_resize_atomic_c_kernel"); 
            int lda = ui_c_get_dim(a).y;
            int ldb = ui_c_get_dim(r).y;
            int size = m*sizeof(T);
            int count = n;
            T* sd = ui_c_current(a)(0,0);
            T* dd = ui_r_updated(r)(0,0);
            do{ memcpy(dd, sd, size); sd += lda; dd += ldb; }while(--count > 0);
            __A_TIME_C_STOP
        }
    };
        
    // {{{ MKL LAPACK kernels

    template<typename T>
    struct svd_atomic : public ambient::kernel_unpinned< svd_atomic<T> > 
    {
        typedef void (svd_atomic::*F)(const maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<T>&, 
                          maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<double>&);

        inline void l(const maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<T>& u, 
                      maquis::types::p_dense_matrix_impl<T>& vt, maquis::types::p_dense_matrix_impl<double>& s)
        {
            this->ctxt_select("1 from ambient as svd_atomic"); //if(!ctxt.involved()) return;
            this->atomic_conditional_assign(ui_l_current(s));
            this->atomic_conditional_assign(ui_l_current(a));
            this->atomic_conditional_assign(ui_l_current(u));
            this->atomic_conditional_assign(ui_l_current(vt));
        }

        inline void c(const maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<T>& u, 
                     maquis::types::p_dense_matrix_impl<T>& vt, maquis::types::p_dense_matrix_impl<double>& s)
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
    struct heev_atomic : public ambient::kernel_unpinned< heev_atomic<T> > 
    {
        typedef void (heev_atomic::*F)(maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<double>&);

        inline void l(maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<double>& w){
            this->ctxt_select("1 from ambient as heev_atomic"); //if(!ctxt.involved()) return;
            this->atomic_conditional_assign(ui_l_current(a));
            this->atomic_conditional_assign(ui_l_current(w));
        }

        inline void c(maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<double>& w){
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
}
#endif
