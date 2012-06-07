#ifndef __MAQUIS_TYPES_KERNELS_ATOMICS_HPP__
#define __MAQUIS_TYPES_KERNELS_ATOMICS_HPP__

extern "C" {
    void dcopy_(const int*, const double*, const int*, double*, const int*);
}


namespace ambient {

    template<typename T>
    struct gemm_general_atomic : public ambient::kernel_dispatch< gemm_general_atomic<T> > 
    { // gs
        typedef void(*F)(const maquis::types::p_dense_matrix_impl<T>&, 
                         const maquis::types::p_dense_matrix_impl<T>&, 
                               maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
            ctxt_select("1 from ambient as gemm_general_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a),0,0);
            atomic_assign(ui_l_current(b),0,0);
            atomic_assign(ui_l_current(c),0,0);
        }
        static inline void c(const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
            __A_TIME("ambient_gemm_general_atomic_c_kernel"); 
            T* bd   = ui_c_current(b)(0,0);
            T* ad   = ui_c_current(a)(0,0);
            T* cd   = ui_c_updated(c)(0,0);
            int m   = ui_c_get_dim(a).y;
            int n   = ui_c_get_dim(b).x;
            int k   = ui_c_get_dim(b).y;
            int lda = ui_c_get_mem_dim(a).y;
            int ldb = ui_c_get_mem_dim(b).y;
            int ldc = ui_c_get_mem_dim(c).y;
            T alpha(1.0); 
            T beta(1.0);
            gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
            __A_TIME_STOP
        }
    };


    template<typename T>
    struct copy_atomic : public ambient::kernel_dispatch< copy_atomic<T> > 
    { // gs
        typedef void(*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& ac, const maquis::types::p_dense_matrix_impl<T>& a){
            ctxt_select("1 from ambient as copy_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a),0,0);
            atomic_assign(ui_l_current(ac),0,0);
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& ac, const maquis::types::p_dense_matrix_impl<T>& a){
            __A_TIME("ambient_copy_atomic_c_kernel"); 
            T* ad  = ui_c_current(a)(0,0);
            T* acd  = ui_c_updated(ac)(0,0);
            __a_copy(acd, ad, ui_c_get_mem_dim(a).square());
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct reshape_l2r_atomic : public ambient::kernel_dispatch< reshape_l2r_atomic<T> > 
    { // gs
        typedef void (*F)(const maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(const maquis::types::p_dense_matrix_impl<T>& left, maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            ctxt_select("1 from ambient as reshape_l2r_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(right),0,0);
            atomic_assign(ui_l_current(left),0,0);
        }

        static inline void c(const maquis::types::p_dense_matrix_impl<T>& left, maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            __A_TIME("ambient_reshape_l2r_atomic_c_kernel"); 
            __a_atomic_refresh(right); // refreshing updated memory
            for(size_t ss = 0; ss < sdim; ++ss){
                __a_memptf_atomic<T, __a_memcpy>(right, dim2(ss*rdim + right_offset, 0), 
                                                 left,  dim2(0, ss*ldim + left_offset), 
                                                 dim2( rdim, ldim ));
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct reshape_r2l_atomic : public ambient::kernel_dispatch< reshape_r2l_atomic<T> > 
    { // gs
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
               const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& left, const maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            ctxt_select("1 from ambient as reshape_l2r_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(left),0,0);
            atomic_assign(ui_l_current(right),0,0);
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& left, const maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            __A_TIME("ambient_reshape_r2l_atomic_c_kernel"); 
            __a_atomic_refresh(left); // refreshing updated memory
            for(size_t ss = 0; ss < sdim; ++ss)
                __a_memptf_atomic<T, __a_memcpy>(left,  dim2(0, ss*ldim + left_offset), 
                                                 right, dim2(ss*rdim + right_offset,0), 
                                                 dim2( rdim, ldim ));
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct lb_tensor_mpo_atomic : public ambient::kernel_dispatch< lb_tensor_mpo_atomic<T> > 
    { // gs
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
               const size_t&, const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            ctxt_select("1 from ambient as rb_tensor_mpo_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(out),0,0);
            atomic_assign(ui_l_current(in),0,0);
            atomic_assign(ui_l_current(alfa),0,0);
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            __A_TIME("ambient_lb_tensor_mpo_atomic_c_kernel"); 
            __a_atomic_refresh(out); // refreshing updated memory
            T* alfad = ui_c_current(alfa)(0,0);
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T alfa_t = alfad[ss1 + ui_c_get_mem_dim(alfa).y*ss2];
                    __a_memptf_atomic<T, __a_memscal>(out, dim2(0, out_offset + ss2*ldim),
                                                      in,  dim2(0, in_offset + ss1*ldim),
                                                      dim2(rdim, ldim), alfa_t);
                }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct rb_tensor_mpo_atomic : public ambient::kernel_dispatch< rb_tensor_mpo_atomic<T> > 
    { // gs
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            ctxt_select("1 from ambient as rb_tensor_mpo_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(out),0,0);
            atomic_assign(ui_l_current(in),0,0);
            atomic_assign(ui_l_current(alfa),0,0);
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            __A_TIME("ambient_rb_tensor_mpo_atomic_c_kernel"); 
            __a_atomic_refresh(out); // refreshing updated memory
            T* alfad = ui_c_current(alfa)(0,0);
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T alfa_t = alfad[ss1 + ui_c_get_mem_dim(alfa).y*ss2];
                    __a_memptf_atomic<T, __a_memscal>(out, dim2(out_offset + ss2*rdim, 0),
                                                      in,  dim2(in_offset + ss1*rdim, 0),
                                                      dim2(rdim, ldim), alfa_t);
                }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct scalar_norm_atomic : public ambient::kernel_dispatch< scalar_norm_atomic<T> > 
    {// gs
        typedef void (*F)(const maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, T*&);

        static inline void l(const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, T*& norm){
            ctxt_select("* from ambient as scalar_norm_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a),0,0);
        }

        static inline void c(const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, T*& norm){
            __A_TIME("ambient_scalar_norm_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            *norm += __a_dot(ad, ad, m*n);;
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct scalar_overlap_atomic : public ambient::kernel_dispatch< scalar_overlap_atomic<T> > 
    { // gs
        typedef void (*F)(const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, T*&);

        static inline void l(const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, const size_t& m, const size_t& n, T*& overlap){
            ctxt_select("* from ambient as scalar_overlap_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a),0,0);
            atomic_assign(ui_l_current(b),0,0);
        }

        static inline void c(const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, const size_t& m, const size_t& n, T*& overlap){
            __A_TIME("ambient_scalar_overlap_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* bd = ui_c_current(b)(0,0);
            *overlap += __a_dot(ad, bd, m*n);;
            __A_TIME_STOP
        }
    };

        
    template<typename T>
    struct add_atomic : public ambient::kernel_dispatch< add_atomic<T> > 
    { // gs
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            ctxt_select("1 from ambient as add_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a),0,0);
            atomic_assign(ui_l_current(b),0,0);
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            __A_TIME("ambient_add_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* bd = ui_c_current(b)(0,0);
            T* ar = ui_c_updated(a)(0,0);
            size_t size = ui_c_get_mem_dim(a).square();
            for(size_t k = 0; k < size; k++)
                ar[k] = ad[k] + bd[k];
            __A_TIME_STOP
        }
    };

        
    template<typename T>
    struct sub_atomic : public ambient::kernel_dispatch< sub_atomic<T> > 
    { // gs
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            ctxt_select("1 from ambient as sub_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a),0,0);
            atomic_assign(ui_l_current(b),0,0);
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            __A_TIME("ambient_sub_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* bd = ui_c_current(b)(0,0);
            T* ar = ui_c_updated(a)(0,0);
            size_t size = ui_c_get_mem_dim(a).square();
            for(size_t k = 0; k < size; k++)
                ar[k] = ad[k] - bd[k];
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct scale_atomic : public ambient::kernel_dispatch< scale_atomic<T> > 
    { // gs
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, const T*&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const T*& t){
            ctxt_select("1 from ambient as scale_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a),0,0);
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const T*& t){
            __A_TIME("ambient_scale_atomic_c_kernel"); 
            T* ad = ui_c_current(a)(0,0);
            T* ar = ui_c_updated(a)(0,0);
            size_t size = m*n;
            for(size_t k=0; k < size; k++) 
                ar[k] = ad[k] * (*t);
            __A_TIME_STOP
        }
    };

    template<typename T>
    struct transpose_out_atomic : public ambient::kernel_dispatch< transpose_out_atomic<T> > 
    { // gs
        typedef void (*F)(const maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&);

        static inline void l(const maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<T>& t, const size_t& m, const size_t& n){
            ctxt_select("1 from ambient as transpose_out_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a),0,0);
            atomic_assign(ui_l_current(t),0,0);
        }

        static inline void c(const maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<T>& t, const size_t& m, const size_t& n){
            __A_TIME("ambient_transpose_out_atomic_c_kernel"); 
            T* od = ui_c_current(a)(0,0);
            T* td = ui_c_updated(t)(0,0);

            size_t mlda = ui_c_get_mem_dim(a).y;
            size_t tlda = ui_c_get_mem_dim(t).y;

            for(size_t j=0; j < n; ++j)
            for(size_t i = 0; i < m; ++i)
            td[j+i*tlda] = od[i+j*mlda];
            
            __A_TIME_STOP
        }
    };
        
    // {{{ MKL LAPACK kernels

    template<typename T>
    struct svd_atomic : public ambient::kernel_dispatch_unpinned< svd_atomic<T> > 
    {
        typedef void (*F)(const maquis::types::p_dense_matrix_impl<T>&, int&, int&, maquis::types::p_dense_matrix_impl<T>&, 
                          maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<double>&);

        static inline void l(const maquis::types::p_dense_matrix_impl<T>& a, int& m, int& n, maquis::types::p_dense_matrix_impl<T>& u, 
                      maquis::types::p_dense_matrix_impl<T>& vt, maquis::types::p_dense_matrix_impl<double>& s)
        {
            ctxt_select("1 from ambient as svd_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(s),0,0);
            atomic_pin(ui_l_current(a),0,0);
            atomic_pin(ui_l_current(u),0,0);
            atomic_pin(ui_l_current(vt),0,0);
        }

        static inline void c(const maquis::types::p_dense_matrix_impl<T>& a, int& m, int& n, maquis::types::p_dense_matrix_impl<T>& u, 
                     maquis::types::p_dense_matrix_impl<T>& vt, maquis::types::p_dense_matrix_impl<double>& s)
        {
            // gs
            __A_TIME("ambient_svd_atomic_c_kernel"); 
        /* Locals */
            int lda = ui_c_get_mem_dim(a).y;
            int ldu = ui_c_get_mem_dim(u).y;
            int ldvt = ui_c_get_mem_dim(vt).y;
            int info, lwork;
            T wkopt;
            T* work;
            double* rwork = new double[5*std::min(lda,ldu)]; // C - useless for double but need for complex 
            T* ad  = ui_c_current(a) (0,0);
            T* ud  = ui_c_updated(u) (0,0);
            T* vtd = ui_c_updated(vt)(0,0);
            double* sd  = ui_c_updated(s) (0,0);
        /* Query and allocate the optimal workspace */
            lwork = -1; // C - Alex, netlib said -1 for the best workspace
            gesvd( "S", "S", &m, &n, ad, &lda, sd, ud, &ldu, vtd, &ldvt, &wkopt, &lwork, rwork, &info );
            lwork = OptimalSize(wkopt);
            work = (T*)malloc( lwork*sizeof(T) );
        /* Compute SVD */
            gesvd( "S", "S", &m, &n, ad, &lda, sd, ud, &ldu, vtd, &ldvt, work, &lwork, rwork, &info );
        /* Check for convergence */
            if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
                exit( 1 );
            }
            free(work);
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct heev_atomic : public ambient::kernel_dispatch_unpinned< heev_atomic<T> > 
    {
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, const size_t&, maquis::types::p_dense_matrix_impl<double>&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, maquis::types::p_dense_matrix_impl<double>& w){
            ctxt_select("1 from ambient as heev_atomic"); //if(!ctxt.involved()) return;
            atomic_pin(ui_l_current(a), 0, 0);
            atomic_pin(ui_l_current(w), 0, 0);
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, maquis::types::p_dense_matrix_impl<double>& w){
            // gs
            __A_TIME("ambient_heev_atomic_c_kernel"); 
            int lda = ui_c_get_mem_dim(a).y;
            int info, lwork = -1;
            double wkopt;
            double* work;
            int am = (int)m; // for mkl (int*)
            double* ad = (double*)__a_solidify<T>(a);
            double* wd = (double*)__a_solidify<T>(w);

            dsyev_("V","U",&am,ad,&lda,wd,&wkopt,&lwork,&info);
            lwork = (int)wkopt;
            work = (double*)malloc( lwork*sizeof(double) );
            dsyev_("V","U",&am,ad,&lda,wd,work,&lwork,&info);
            if( info > 0 ) {
                printf( "The algorithm computing SYEV failed to converge.\n" );
                exit( 1 );
            }
            // First we reverse the eigenvalues, to be in agreement with the serial version ! 
            // The matrix is solidified, so we do not care on the workgroup representation
            double tempdbl;
            for (int i=0; i< static_cast<int>(m/2); i++){ 
                tempdbl = wd[i];
                wd[i] = wd[m-i-1];
                wd[m-i-1] = tempdbl;
            } 
            // Second we reverse the eigenvectors
            double* tempcol = new double[lda]; 
            for (int i=0; i< static_cast<int>(m/2); ++i){ 
                memmove((void*)tempcol,(void*)&ad[i*lda],lda*sizeof(double));
                memmove((void*)&ad[i*lda],(void*)&ad[(m-1-i)*lda],lda*sizeof(double));
                memmove((void*)&ad[(m-1-i)*lda],(void*)tempcol,lda*sizeof(double));
            }
            delete[] tempcol; 
         
            __a_disperse<T>(ad, a);
            __a_disperse<T>(wd, w);
            free(work);
            __A_TIME_STOP
        }
    };

    // }}}
}
#endif
