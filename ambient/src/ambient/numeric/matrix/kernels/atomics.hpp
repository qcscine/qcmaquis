#ifndef __AMBIENT_NUMERIC_MATRIX_KERNELS_ATOMICS_HPP__
#define __AMBIENT_NUMERIC_MATRIX_KERNELS_ATOMICS_HPP__

extern "C" {
    void dgemm_(const char*,const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, double*, const int*);
}

namespace ambient { namespace numeric { namespace kernels {

    using ambient::numeric::matrix_impl;

    template<class ViewA, class ViewB, typename T>
    struct gemm_general_atomic : public kernel_atomic< gemm_general_atomic<ViewA, ViewB, T> > 
    { // gs
        typedef void(gemm_general_atomic::*F)(const matrix_impl<T>&, const matrix_impl<T>&, weak_matrix_impl<T>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<T>& b, weak_matrix_impl<T>& c){
            this->ctxt_select("1 from ambient as gemm_general_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b));
            this->assign(ui_l_current(c));
        }
        inline void c(const matrix_impl<double>& a, const matrix_impl<double>& b, weak_matrix_impl<double>& c){
            __A_TIME_C("ambient_gemm_general_atomic_c_kernel"); 
            double* ad = ui_c_current(a);
            double* bd = ui_c_current(b);
            double* cd = ui_w_updated(c);
            int m = ViewA::rows(a);
            int k = ViewA::cols(a);
            int n = ViewB::cols(b);
            int lda = ui_c_get_dim(a).y;
            int ldb = ui_c_get_dim(b).y;
            int ldc = ui_c_get_dim(c).y;
            static const double alpha(1.0); 
            static const double beta(0.0);
            dgemm_(ViewA::code(), ViewB::code(), &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
            __A_TIME_C_STOP
        }
        inline void c(const matrix_impl<std::complex<double> >& a, const matrix_impl<std::complex<double> >& b, weak_matrix_impl<std::complex<double> >& c){
            __A_TIME_C("ambient_gemm_general_atomic_c_kernel"); 
            T* ad   = ui_c_current(a);
            T* bd   = ui_c_current(b);
            T* cd   = ui_w_updated(c);
            int m   = ui_c_get_dim(a).y;
            int n   = ui_c_get_dim(b).x;
            int k   = ui_c_get_dim(b).y;
            T alpha(1.0); 
            T beta(0.0);
            gemm("N","N", &m, &n, &k, &alpha, ad, &m, bd, &k, &beta, cd, &m);
            __A_TIME_C_STOP
        }
    };
        
    template<class ViewB, typename T, typename D>
    struct gemm_diagonal_lhs : public kernel_atomic< gemm_diagonal_lhs<ViewB,T,D> > 
    {
        typedef void (gemm_diagonal_lhs::*F)(const matrix_impl<D>&, const matrix_impl<T>&, weak_matrix_impl<T>&);

        inline void l(const matrix_impl<D>& a_diag, const matrix_impl<T>& b, weak_matrix_impl<T>& c){
            this->ctxt_select("1 from ambient as gemm_diagonal_lhs"); //if(!ctxt.involved()) return;
            this->assign(ui_l_current(a_diag));
            this->pin(ui_l_current(b));
            this->assign(ui_l_current(c));
        }

        inline void c(const matrix_impl<D>& a_diag, const matrix_impl<T>& b, weak_matrix_impl<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_diagonal_lhs_c_kernel"); 
            int sizey = ui_c_get_dim(a_diag).y;
            int size = ui_c_get_dim(b).x;
            int ONE  = 1;
            D* bd = ui_c_current(b);
            D* cd = ui_p_updated(c);
            D* alpha = ui_c_current(a_diag);
        
            for(int k = 0 ; k < sizey; k++){
        	     axpy(&size, &alpha[k], &bd[k], &sizey, &cd[k], &sizey);
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T, typename D>
    struct gemm_diagonal_lhs<transpose_view<matrix<T> >,T,D> : public kernel_atomic< gemm_diagonal_lhs<transpose_view<matrix<T> >,T,D> > 
    {
        typedef void (gemm_diagonal_lhs::*F)(const matrix_impl<D>&, const matrix_impl<T>&, weak_matrix_impl<T>&);

        inline void l(const matrix_impl<D>& a_diag, const matrix_impl<T>& b, weak_matrix_impl<T>& c){
            this->ctxt_select("1 from ambient as gemm_diagonal_lhs"); //if(!ctxt.involved()) return;
            this->assign(ui_l_current(a_diag));
            this->pin(ui_l_current(b));
            this->assign(ui_l_current(c));
        }

        inline void c(const matrix_impl<D>& a_diag, const matrix_impl<T>& b, weak_matrix_impl<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_diagonal_lhs_c_kernel"); 
            printf("Special DIAGONAL!\n");
            size_t sizex = ui_c_get_dim(b).x;
            int size  = ui_c_get_dim(a_diag).y;
            int ONE  = 1;
            D* bd = ui_c_current(b);
            D* cd = ui_p_updated(c);
            D* alpha = ui_c_current(a_diag);
        
            for(int k = 0 ; k < sizex; k++){
        	     axpy(&size, &alpha[k], &bd[k*size], &ONE, &cd[k], &size);
            }
            __A_TIME_C_STOP
        }
    };
        
    template<class ViewA, typename T, typename D>
    struct gemm_diagonal_rhs : public kernel_atomic< gemm_diagonal_rhs<ViewA,T,D> > 
    {
        typedef void (gemm_diagonal_rhs::*F)(const matrix_impl<T>&, const matrix_impl<D>&, weak_matrix_impl<T>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<D>& b_diag, weak_matrix_impl<T>& c){
            this->ctxt_select("1 from ambient as gemm_diagonal_rhs"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b_diag));
            this->assign(ui_l_current(c));
        }

        inline void c(const matrix_impl<T>& a, const matrix_impl<D>& b_diag, weak_matrix_impl<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_diagonal_rhs_c_kernel"); 
            size_t sizex = ui_c_get_dim(b_diag).y;
            int size = ui_c_get_dim(a).y; // for the case of complex
            int ONE = 1;
            D* ad = ui_c_current(a);
            D* cd = ui_p_updated(c);
        	D* alpha = ui_c_current(b_diag);
        
            for(int k = 0 ; k < sizex; k++){
        	    axpy(&size, &alpha[k], &ad[k*size], &ONE, &cd[k*size], &ONE);
            }
            __A_TIME_C_STOP
        }
    };

    template<typename T, typename D>
    struct gemm_diagonal_rhs<transpose_view<matrix<T> >,T,D> : public kernel_atomic< gemm_diagonal_rhs<transpose_view<matrix<T> >,T,D> > 
    {
        typedef void (gemm_diagonal_rhs::*F)(const matrix_impl<T>&, const matrix_impl<D>&, weak_matrix_impl<T>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<D>& b_diag, weak_matrix_impl<T>& c){
            this->ctxt_select("1 from ambient as gemm_diagonal_rhs"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(b_diag));
            this->assign(ui_l_current(c));
        }

        inline void c(const matrix_impl<T>& a, const matrix_impl<D>& b_diag, weak_matrix_impl<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_diagonal_rhs_c_kernel"); 
            printf("Special DIAGONAL!\n");
            int sizey = ui_c_get_dim(b_diag).y;
            int size = ui_c_get_dim(a).x;
            int ONE = 1;
            D* ad = ui_c_current(a);
            D* cd = ui_p_updated(c);
        	D* alpha = ui_c_current(b_diag);
        
            for(int k = 0 ; k < sizey; k++){
        	    axpy(&size, &alpha[k], &ad[k], &sizey, &cd[k*size], &ONE);
            }
            __A_TIME_C_STOP
        }
    };


    template<typename T>
    struct copy_atomic : public kernel_atomic< copy_atomic<T> > 
    { // gs
        typedef void(copy_atomic::*F)(weak_matrix_impl<T>&, const matrix_impl<T>&);

        inline void l(weak_matrix_impl<T>& ac, const matrix_impl<T>& a){
            this->ctxt_select("1 from ambient as copy_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(ac));
        }

        inline void c(weak_matrix_impl<T>& ac, const matrix_impl<T>& a){
            __A_TIME_C("ambient_copy_atomic_c_kernel"); 
            T* ad  = ui_c_current(a);
            T* acd  = ui_w_updated(ac);
            memcpy(acd, ad, ui_c_get_mem_size(a));
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct op_kron_atomic : public kernel_atomic< op_kron_atomic<T> > 
    { // gs - 2su
        typedef void (op_kron_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&, const matrix_impl<T>&,
                                          const size_t&, const size_t&, 
                                          const size_t&, const size_t&,
                                          const size_t&, const size_t&);

        inline void l(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_y_offset, const size_t& out_x_offset, 
                      const size_t& ldim1, const size_t& ldim2, 
                      const size_t& rdim1, const size_t& rdim2)
        {
            this->ctxt_select("1 from ambient as op_kron_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(out));
            this->assign(ui_l_current(in));
            this->assign(ui_l_current(alfa));
        }

        inline void c(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_y_offset, const size_t& out_x_offset, 
                      const size_t& ldim1, const size_t& ldim2, 
                      const size_t& rdim1, const size_t& rdim2)
        {
            __A_TIME_C("ambient_op_kron_atomic_c_kernel"); 
            T* alfad = ui_c_current(alfa);
            for(size_t l1 = 0; l1 < ldim1; ++l1)
            for(size_t r1 = 0; r1 < rdim1; ++r1)
            {
                T alfa_t = alfad[l1 + r1*ui_c_get_dim(alfa).y];
                __a_memptf_atomic_r<T, __a_memscal>(out, dim2(out_x_offset + r1*rdim2, out_y_offset + l1*ldim2),
                                                    in,  dim2(0, 0), dim2(rdim2, ldim2), alfa_t);
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct reshape_l2b_atomic : public kernel_atomic< reshape_l2b_atomic<T> > 
    { // gs - 2su
        typedef void (reshape_l2b_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&,
                                              const size_t&, const size_t&, 
                                              const size_t&, const size_t&, 
                                              const size_t&, const size_t&,
                                              const size_t&, const size_t&);

        inline void l(matrix_impl<T>& out, const matrix_impl<T>& in,
                      const size_t& in_left_offset, const size_t& in_phys_offset, 
                      const size_t& out_left_offset, const size_t& out_x_offset,
                      const size_t& sdim1, const size_t& sdim2, 
                      const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as reshape_l2b_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(out));
            this->assign(ui_l_current(in));
        }

        inline void c(matrix_impl<T>& out, const matrix_impl<T>& in,
                      const size_t& in_left_offset, const size_t& in_phys_offset, 
                      const size_t& out_left_offset, const size_t& out_x_offset,
                      const size_t& sdim1, const size_t& sdim2, 
                      const size_t& ldim, const size_t& rdim)
        {
            __A_TIME_C("ambient_reshape_l2b_atomic_c_kernel"); 

            size_t in_y_offset  = in_left_offset + ldim*in_phys_offset;
            size_t out_y_offset = out_left_offset;

            __a_atomic_refresh(out);
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1){
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    __a_memptf_atomic_r<T, __a_memcpy>(out, dim2(out_x_offset + rdim*ss2, out_y_offset), 
                                                       in,  dim2(0, in_y_offset), 
                                                       dim2( rdim, ldim ));
                    in_y_offset += ldim;
                }
                out_y_offset += ldim;
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct reshape_b2l_atomic : public kernel_atomic< reshape_b2l_atomic<T> > 
    { // gs - 2su
        typedef void (reshape_b2l_atomic::*F)(matrix_impl<T>&, const matrix_impl<T>&,
                                              const size_t&, const size_t&, 
                                              const size_t&, const size_t&, 
                                              const size_t&, const size_t&,
                                              const size_t&, const size_t&);

        inline void l(matrix_impl<T>& out, const matrix_impl<T>& in,
                      const size_t& in_left_offset, const size_t& in_x_offset, 
                      const size_t& out_left_offset, const size_t& out_phys_offset,
                      const size_t& sdim1, const size_t& sdim2, 
                      const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as reshape_b2l_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(out));
            this->assign(ui_l_current(in));
        }

        inline void c(matrix_impl<T>& out, const matrix_impl<T>& in,
                      const size_t& in_left_offset, const size_t& in_x_offset, 
                      const size_t& out_left_offset, const size_t& out_phys_offset,
                      const size_t& sdim1, const size_t& sdim2, 
                      const size_t& ldim, const size_t& rdim)
        {
            __A_TIME_C("ambient_reshape_b2l_atomic_c_kernel"); 

            size_t in_y_offset  = in_left_offset;
            size_t out_y_offset = out_left_offset + out_phys_offset*ldim;

            __a_atomic_refresh(out);
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1){
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2)
                {
                    __a_memptf_atomic_r<T, __a_memcpy>(out, dim2(0, out_y_offset), 
                                                       in,  dim2(in_x_offset + rdim*ss2, in_y_offset), 
                                                       dim2( rdim, ldim ));
                    out_y_offset += ldim;
                }
                in_y_offset += ldim;
            }
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
            T* alfad = ui_c_current(alfa);
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T alfa_t = alfad[ss1 + ss2*ui_c_get_dim(alfa).y];
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
            T* alfad = ui_c_current(alfa);
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
    struct trace_atomic : public kernel_atomic< trace_atomic<T> > 
    {
        typedef void (trace_atomic::*F)(const matrix_impl<T>&, future<T>&);

        inline void l(const matrix_impl<T>& a, future<T>& trace){
            this->ctxt_select("* from ambient as trace_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(const matrix_impl<T>& a, future<T>& trace){
            // gs
            __A_TIME_C("ambient_trace_atomic_c_kernel"); 
            size_t m = ui_c_get_dim(a).y;
            size_t n = ui_c_get_dim(a).x;
            T* ad = ui_c_current(a);
        
            size_t sizex = std::min(n,m);
            for(size_t jj = 0; jj < sizex; jj++){
                trace.get_value() += ad[jj + jj*m];
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
            T* ad = ui_c_current(a);
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
            T* ad = ui_c_current(a);
            T* bd = ui_c_current(b);
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
            T* ad = ui_c_current(a);
            T* bd = ui_c_current(b);
            T* ar = ui_r_updated(a);
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
            T* bd = ui_c_current(b);
            T* ar = ui_r_updated(a);
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
            T* ar = ui_r_updated(a);
            int size = ui_c_get_dim(a).square();
            static const int ONE = 1;
            dscal_( &size, &t.get_value(), ar, &ONE );
            __A_TIME_C_STOP
        }

        inline void c(matrix_impl<std::complex<double> >& a, const future< std::complex<double> >& t){
            __A_TIME_C("ambient_scale_atomic_c_kernel"); 
            T* ad = ui_c_current(a);
            T* ar = ui_w_updated(a);
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
            T* ar = ui_r_updated(a);
            int size = ui_c_get_dim(a).square();
            static const int ONE = 1;
            double factor = 1. / t.get_value();
            dscal_( &size, &factor, ar, &ONE );
            __A_TIME_C_STOP
        }

        inline void c(matrix_impl<std::complex<double> >& a, const future< std::complex<double> >& t){
            __A_TIME_C("ambient_scale_inverse_atomic_c_kernel"); 
            T* ad = ui_c_current(a);
            T* ar = ui_w_updated(a);
            int size = ui_c_get_dim(a).square();
            for(int k=0; k < size; k++) 
                ar[k] = ad[k] / t.get_value();
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct transpose_out_atomic : public kernel_atomic< transpose_out_atomic<T> > 
    { // gs
        typedef void (transpose_out_atomic::*F)(const matrix_impl<T>&, weak_matrix_impl<T>&);

        inline void l(const matrix_impl<T>& a, weak_matrix_impl<T>& t){
            this->ctxt_select("1 from ambient as transpose_out_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(t));
        }

        inline void c(const matrix_impl<T>& a, weak_matrix_impl<T>& t){
            __A_TIME_C("ambient_transpose_out_atomic_c_kernel"); 
            T* od = ui_c_current(a);
            T* td = ui_w_updated(t);
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
        typedef void (resize_atomic::*F)(weak_matrix_impl<T>&, const matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(weak_matrix_impl<T>& r, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            this->ctxt_select("1 from ambient as resize_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(r));
        }

        inline void c(weak_matrix_impl<T>& r, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            __A_TIME_C("ambient_resize_atomic_c_kernel"); 
            __a_memptf_atomic_r<T, __a_memcpy>(r, dim2(0,0), a, dim2(0,0), dim2(n, m)); 
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct init_identity_atomic : public kernel_atomic< init_identity_atomic<T> > 
    {
        typedef void (init_identity_atomic::*F)(weak_matrix_impl<T>&);

        inline void l(weak_matrix_impl<T>& a){
            this->ctxt_select("1 from ambient as init_identity_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(weak_matrix_impl<T>& a){
            __A_TIME_C("ambient_init_identity_atomic_c_kernel"); 
            size_t n = ui_c_get_dim(a).x;
            size_t m = ui_c_get_dim(a).y;
            T* ad = ui_r_updated(a);

            size_t sizex = std::min(m,n); // respecting borders
            for(size_t jj = 0; jj < sizex; ++jj) ad[jj + m*jj] = 1.;
            __A_TIME_C_STOP
        }
    };
       
    template<typename T>
    struct init_random_atomic : public kernel_atomic< init_random_atomic<T> > 
    {
        typedef void (init_random_atomic::*F)(weak_matrix_impl<T>&);
     
        template<typename T> inline void randomize(T* ad){ *ad = drand48(); }
        template<typename T> inline void randomize(std::complex<T>* ad){
            ad->real(drand48());
            ad->imag(drand48());
        }

        inline void l(weak_matrix_impl<T>& a){
            this->ctxt_select("1 from ambient as init_random_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }
        
        inline void c(weak_matrix_impl<T>& a){
            __A_TIME_C("ambient_init_random_atomic_c_kernel"); 
            size_t m = ui_c_get_dim(a).y;
            size_t n = ui_c_get_dim(a).x;
            T* ad = ui_w_updated(a);
          
            for(size_t jj = 0; jj < n; jj++){
                for(size_t ii = 0; ii < m; ii++)
                    randomize((ad+(jj*m+ii)));
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct init_value_atomic : public kernel_atomic< init_value_atomic<T> > 
    {
        typedef void (init_value_atomic::*F)(weak_matrix_impl<T>&, const T&);

        inline void l(weak_matrix_impl<T>& a, const T& value){
            this->ctxt_select("1 from ambient as init_value_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(weak_matrix_impl<T>& a, const T& value){
            __A_TIME_C("ambient_init_value_atomic_c_kernel"); 
            T* ad = ui_w_updated(a);
            size_t m = ui_c_get_dim(a).y;
            size_t n = ui_c_get_dim(a).x;
        
            for(size_t j=0; j < n; ++j){
                for(size_t i = 0; i < m; ++i){
                    ad[i+j*m] = value; // not a memset due to complex
                }
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct round_square_atomic : public kernel_atomic< round_square_atomic<T> > 
    {
        typedef void (round_square_atomic::*F)(const matrix_impl<T>&, std::vector<T>*&);

        inline void l(const matrix_impl<T>& a, std::vector<T>*& ac){
            this->ctxt_select("* from ambient as round_square_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(const matrix_impl<T>& a, std::vector<T>*& ac){
            // gs
            __A_TIME_C("ambient_round_square_atomic_c_kernel"); 
            T* ad = ui_c_current(a);
            size_t sizey = ui_c_get_dim(a).y;
            for(int i=0; i < sizey; i++){
                double v = std::abs(ad[i]);
                if(v > 1e-10) ac->push_back(v*v);
            }
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct cast_to_vector_atomic : public kernel_atomic< cast_to_vector_atomic<T> > 
    {
        typedef void (cast_to_vector_atomic::*F)(std::vector<T>*&, const matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(std::vector<T>*& ac, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            this->ctxt_select("* from ambient as cast_to_vector_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(std::vector<T>*& ac, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            // gs
            __A_TIME_C("ambient_cast_to_vector_atomic_c_kernel"); 
            T* ad = ui_c_current(a);
            for(int j=0; j < n; ++j) memcpy((void*)&(*ac)[j*m],(void*)&ad[j*m], m*sizeof(T));  
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct cast_from_vector_atomic : public kernel_atomic< cast_from_vector_atomic<T> > 
    {
        typedef void (cast_from_vector_atomic::*F)(const std::vector<T>*&, matrix_impl<T>&, const size_t&, const size_t&, const size_t&);

        inline void l(const std::vector<T>*& ac, matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& lda){
            this->ctxt_select("1 from ambient as cast_from_vector"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
        }

        inline void c(const std::vector<T>*& ac, matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& lda){
            __A_TIME_C("ambient_cast_from_vector_atomic_c_kernel"); 
            T* ad = ui_w_updated(a);
            for(int j=0; j < n; ++j) memcpy((void*)&ad[j*m],(void*)&(*ac)[j*lda], m*sizeof(T));
            __A_TIME_C_STOP 
        }
    };

    template<typename T>
    struct validation_atomic : public kernel_atomic< validation_atomic<T> > 
    {
        typedef void (validation_atomic::*F)(const matrix_impl<T>&, const matrix_impl<T>&, future<int>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<T>& b, future<int>& ret){
            this->ctxt_select("1 from ambient as validation_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a)); 
            this->assign(ui_l_current(b)); 
        }
        
        inline void c(const matrix_impl<T>& a, const matrix_impl<T>& b, future<int>& ret){ // see paper for Reference Dongara 
            T* ad = ui_c_current(a); 
            T* bd = ui_c_current(b); 
            double res; 
            double epsilon = std::numeric_limits<double>::epsilon();
            size_t position_xy; 
       
            for(size_t ii=0; ii < ui_c_get_dim(a).y; ++ii){
                for(size_t jj=0; jj < ui_c_get_dim(a).x; ++jj){
                    if(jj < std::min(ui_c_get_dim(a).x,ui_c_get_dim(b).x) && 
                       ii < std::min(ui_c_get_dim(a).y,ui_c_get_dim(b).y)  )
                    {
                        position_xy = jj*ui_c_get_dim(a).y+ii;
                        res = (norm(ad[position_xy])-norm(bd[position_xy]))/fabs(epsilon*norm(bd[position_xy])); // to do : rotation pb  with complex to change
                        if(res > 256){ // 16 is recommended by Dongara, 256 because lapack gives != runs after runs
                            std::cout << ii << " " << jj << " : " << ad[position_xy] << " " << bd[position_xy] << std::endl;
                            ret.get_value() = 0; // test failed return 0 (bool false)
                        }
                    }
                }
            }
        }
    };

    // {{{ MKL LAPACK kernels

    template<typename T>
    struct svd_atomic : public kernel_atomic< svd_atomic<T> > 
    {
        typedef void (svd_atomic::*F)(const matrix_impl<T>&, weak_matrix_impl<T>&, weak_matrix_impl<T>&, weak_matrix_impl<double>&);

        inline void l(const matrix_impl<T>& a, weak_matrix_impl<T>& u, weak_matrix_impl<T>& vt, weak_matrix_impl<double>& s)
        {
            this->ctxt_select("1 from ambient as svd_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(s));
            this->assign(ui_l_current(u));
            this->assign(ui_l_current(vt));
        }

        inline void c(const matrix_impl<T>& a, weak_matrix_impl<T>& u, weak_matrix_impl<T>& vt, weak_matrix_impl<double>& s)
        { // gs
            __A_TIME_C("ambient_svd_atomic_c_kernel"); 
            int m = ui_c_get_dim(a).y;
            int n = ui_c_get_dim(a).x;
            int k = std::min(m,n);
            int info;
            int lwork = -1; // C - Alex, netlib said -1 for the best workspace
            T wkopt;
            T* ad  = ui_c_current(a);
            T* ud  = ui_r_updated(u);
            T* vtd = ui_r_updated(vt);
            double* sd  = ui_r_updated(s);
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
    struct heev_atomic : public kernel_atomic< heev_atomic<T> > 
    {
        typedef void (heev_atomic::*F)(matrix_impl<T>&, weak_matrix_impl<double>&);

        inline void l(matrix_impl<T>& a, weak_matrix_impl<double>& w){
            this->ctxt_select("1 from ambient as heev_atomic"); //if(!ctxt.involved()) return;
            this->pin(ui_l_current(a));
            this->assign(ui_l_current(w));
        }

        inline void c(matrix_impl<T>& a, weak_matrix_impl<double>& w){
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
