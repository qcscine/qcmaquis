#ifndef __MAQUIS_TYPES_KERNELS_HPP__
#define __MAQUIS_TYPES_KERNELS_HPP__

namespace ambient {

    template<typename T>
    struct gemm_inplace : public ambient::kernel_dispatch< gemm_inplace<T> > 
    {
        typedef void(*F)(pinned maquis::types::p_dense_matrix_impl<T>&, 
                          const maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            ctxt_select("1 from ambient as gemm_inplace"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
            block_2d_cycle_assign(ui_l_current(b));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            int m   = ui_c_get_mem_dim(a).y;
            int n   = ui_c_get_mem_dim(b).x;
            int k   = ui_c_get_mem_dim(b).y;
            int lda = m;
            int ldb = k;
            int ldc = m;
            T alpha(1.0); 
            T beta(1.0);
        // a(x,y) => a(x,z) x b(y,x)  where z : [m,0)
        // current block of matrix a:
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
        // taking (y,x) of b:
            if(ui_c_get_grid_dim(b).y > x) while(y < ui_c_get_grid_dim(b).x){
                T* bd = ui_c_current(b)(y,x); // remote
        // multiplying with column of a:
                std::list<int> L;
                for(int z = 0; z < ui_c_get_grid_dim(a).y; z++) L.push_back(z);
                while(!L.empty()){
                    std::list<int>::iterator zy = L.begin();
                    while(zy != L.end()){
                        if(!ui_c_updated(a)(y,*zy).trylock()){ zy++; continue; }
                        T* ad = ui_c_current(a)(x,*zy);
                        T* cd = ui_c_updated(a)(y,*zy); // a(x,z) x b(y,x) => c(y,z)
                        gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
                        ui_c_updated(a)(y,*zy).unlock();
                        L.erase(zy++);
                    }
                }
                y += ui_c_get_grid_dim(a).y;
            }
        }
    };

    template<typename T>
    struct gemm_general : public ambient::kernel_dispatch< gemm_general<T> > 
    {
        typedef void(*F)(pinned const maquis::types::p_dense_matrix_impl<T>&, 
                                const maquis::types::p_dense_matrix_impl<T>&, 
                                      maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
            ctxt_select("1 from ambient as gemm"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
            block_2d_cycle_assign(ui_l_current(b));
            block_2d_cycle_assign(ui_l_current(c));
        }

        static inline void c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
            // gs
            __A_TIME("ambient_gemm_c_kernel");
            //if(ui_c_get_grid_dim(a) == 1 && ui_c_get_grid_dim(b) == 1) // use gemm atomic for this case
            if(ui_c_get_mem_dim(a) != ui_c_get_mem_dim(b) || ui_c_get_mem_dim(a) != ui_c_get_mem_dim(c)){
                T* ad = (T*)__a_solidify(a);
                T* bd = (T*)__a_solidify(b);
                T* cd = (T*)__a_solidify(c);
                int m   = ui_c_get_dim(a).y;
                int n   = ui_c_get_dim(b).x;
                int k   = ui_c_get_dim(b).y;
                int lda = ui_c_get_mem_dim(a).y*ui_c_get_grid_dim(a).y;
                int ldb = ui_c_get_mem_dim(b).y*ui_c_get_grid_dim(b).y;
                int ldc = ui_c_get_mem_dim(c).y*ui_c_get_grid_dim(c).y;
                T alpha(1.0); 
                T beta(1.0);
                gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
                __a_disperse(cd, c);
                __A_TIME_STOP
                return;
            }
            int m   = ui_c_get_mem_dim(a).y;
            int n   = ui_c_get_mem_dim(b).x;
            int k   = ui_c_get_mem_dim(b).y;
            int lda = m;
            int ldb = k;
            int ldc = m;
            T alpha(1.0); 
            T beta(1.0);
        // a(x,y) => a(x,z) x b(y,x)  where z : [m,0)
        // current block of matrix a:
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
        
            dim2 a_grid_dim = ui_c_get_grid_dim(a);
            dim2 b_grid_dim = ui_c_get_grid_dim(b);
        // taking (y,x) of b:
            if(b_grid_dim.y > x) while(y < b_grid_dim.x){
                T* bd = ui_c_current(b)(y,x); // remote
        // multiplying with column of a:
                std::list<int> L;
                for(int z = 0; z < a_grid_dim.y; z++) L.push_back(z);
                while(!L.empty()){
                    std::list<int>::iterator zy = L.begin();
                    while(zy != L.end()){
                        if(!ui_c_updated(c)(y,*zy).trylock()){ zy++; continue; }
                        T* ad = ui_c_current(a)(x,*zy);
                        T* cd = ui_c_updated(c)(y,*zy); // a(x,z) x b(y,x) => c(y,z)
                        gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
                        ui_c_updated(c)(y,*zy).unlock();
                        L.erase(zy++);
                    }
                }
                y += a_grid_dim.y;
            }
            __A_TIME_STOP
        }
    };


    template<typename T>
    struct copy : public ambient::kernel_dispatch< copy<T> > 
    {
        typedef void(*F)(maquis::types::p_dense_matrix_impl<T>&, pinned const maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a){
            ctxt_select("1 from ambient as copy"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
            block_2d_cycle_assign(ui_l_current(ac));
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a){
            // gs
            __A_TIME("ambient_copy_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            T* a_elements  = ui_c_current(a)(x,y);
            T* ac_elements = ui_c_updated(ac)(x,y);
            memcpy(ac_elements, a_elements, sizeof(T)*ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y);
            __A_TIME_STOP
        }
    };

    template<typename T>
    struct remove_rows : public ambient::kernel_dispatch< remove_rows<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& i_mark, const size_t& k){
            ctxt_select("1 from ambient as remove_rows"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }
        
        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& i_mark, const size_t& k){
            size_t numcols = ui_c_get_dim(a).x;
            size_t numrows = ui_c_get_dim(a).y;
            __a_memptf_reverse<T, __a_memcpy>(a, dim2(0,0), a, dim2(0,0), dim2(numcols, i_mark));
            __a_memptf_reverse<T, __a_memcpy>(a, dim2(0,i_mark), a, dim2(0,k+i_mark), dim2(numcols,numrows-k-i_mark));
        }
    };
        
    template<typename T>
    struct remove_cols : public ambient::kernel_dispatch< remove_cols<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& j_mark, const size_t& k){
            ctxt_select("1 from ambient as remove_cols"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& j_mark, const size_t& k){
            size_t numcols = ui_c_get_dim(a).x;
            size_t numrows = ui_c_get_dim(a).y;
            __a_memptf_reverse<T, __a_memcpy>(a, dim2(0,0), a, dim2(0,0), dim2(j_mark, numrows));
            __a_memptf_reverse<T, __a_memcpy>(a, dim2(j_mark,0), a, dim2(k+j_mark,0), dim2(numcols-k-j_mark,numrows));
        }
    };
        
        
    template<typename T>
    struct resize : public ambient::kernel_dispatch< resize<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& om, const size_t& on){
            ctxt_select("1 from ambient as resize"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& om, const size_t& on){
            __a_memptf_reverse<T, __a_memcpy>(a, dim2(0,0), a, dim2(0,0), dim2(std::min(n,on), std::min(m,om)));
        }
    };
        
    template<typename T>
    struct sqrt_diagonal : public ambient::kernel_dispatch< sqrt_diagonal<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a){
            ctxt_select("1 from ambient as sqrt_diagonal"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a){
            T* ad = ui_c_current(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            T* sd = ui_c_updated(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            size_t size = ui_c_get_mem_dim(a).y;
            for(int i=0; i < size; i++)
                sd[i] = sqrt(ad[i]);
        }
    };
        
    template<typename T>
    struct exp_diagonal : public ambient::kernel_dispatch< exp_diagonal<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const T&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const T& alfa){
            ctxt_select("1 from ambient as exp_diagonal"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const T& alfa){
            T* ad = ui_c_current(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            T* sd = ui_c_updated(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            size_t size = ui_c_get_mem_dim(a).y;
            for(int i=0; i < size; i++)
                sd[i] = exp(ad[i]*alfa);
        }
    };
        
    template<typename T>
    struct exp_diagonal_rc : public ambient::kernel_dispatch< exp_diagonal_rc<T> > 
    {
        typedef void (*F)(maquis::types::p_dense_matrix_impl< std::complex<T> >&, pinned const maquis::types::p_dense_matrix_impl<T>&, const std::complex<T>&);

        static inline void l(maquis::types::p_dense_matrix_impl< std::complex<T> >& e, pinned const maquis::types::p_dense_matrix_impl<T>& a, const std::complex<T>& alfa){
            ctxt_select("1 from ambient as exp_diagonal"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(e));
            block_2d_cycle_assign(ui_l_current(a));
        }

        static inline void c(maquis::types::p_dense_matrix_impl< std::complex<T> >& e, pinned const maquis::types::p_dense_matrix_impl<T>& a, const std::complex<T>& alfa){
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            std::complex<T>* ed = ui_c_updated(e)(x,y);
            T* ad = ui_c_current(a)(x,y);
            size_t size = ui_c_get_mem_dim(e).y;
            for(int i=0; i < size; i++)
                ed[i] = exp(ad[i]*alfa);
        }
    };
        
    template<typename T>
    struct push_back_sqr_gt : public ambient::kernel_dispatch< push_back_sqr_gt<T> > 
    {
        typedef void (*F)(pinned const maquis::types::p_dense_matrix_impl<T>&, std::vector<T>*&);

        static inline void l(pinned const maquis::types::p_dense_matrix_impl<T>& a, std::vector<T>*& ac){
            ctxt_select("* from ambient as push_back_sqr_gt"); //if(!ctxt.involved()) return;
            block_outright_pin(ui_l_current(a));
        }

        static inline void c(pinned const maquis::types::p_dense_matrix_impl<T>& a, std::vector<T>*& ac){
            // gs
            __A_TIME("ambient_push_back_sqr_gt_c_kernel"); 
            T* ad = ui_c_current(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            size_t sizey = __a_get_limit_y(a);
            for(int i=0; i < sizey; i++){
                double v = std::abs(ad[i]);
                if(v > 1e-10) ac->push_back(v*v);
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct cast_to_dense : public ambient::kernel_dispatch< cast_to_dense<T> > 
    {
        typedef void (*F)(std::vector<T>*&, pinned const maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&);

        static inline void l(std::vector<T>*& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
            ctxt_select("* from ambient as cast_to_dense"); //if(!ctxt.involved()) return;
            block_outright_pin(ui_l_current(a));
        }

        static inline void c(std::vector<T>*& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
            // gs
            __A_TIME("ambient_cast_to_dense_c_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            int xx = ui_c_get_mem_dim(a).x*x; 
            int yy = ui_c_get_mem_dim(a).y*y; // conversion cartersia coordinates dense / p_dense
            size_t offset;
            T* ad = ui_c_current(a)(x,y);
           
            size_t sizex = __a_get_limit_x(a, n);
            size_t sizey = __a_get_limit_y(a, m);
            size_t lda = ui_c_get_mem_dim(a).y;
        
            for(int j=0; j < sizex; ++j){
                offset = yy + (xx+j)*m;
                memcpy((void*)&(*ac)[offset],(void*)&ad[j*lda], sizey*sizeof(T));  
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct cast_to_p_dense : public ambient::kernel_dispatch< cast_to_p_dense<T> > 
    {
        typedef void (*F)(const std::vector<T>*&, pinned maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, const size_t&);

        static inline void l(const std::vector<T>*& ac, pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& lda){
            ctxt_select("1 from ambient as cast_to_p_dense"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }

        static inline void c(const std::vector<T>*& ac, pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& lda){
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            size_t xx = ui_c_get_mem_dim(a).x*x;
            size_t yy = ui_c_get_mem_dim(a).y*y;
            size_t sizex = __a_get_limit_x(a, n);
            size_t sizey = __a_get_limit_y(a, m);

            T* ad = ui_c_updated(a)(x,y);
            for(int j=0; j < sizex; ++j)
                memcpy((void*)&ad[j*ui_c_get_mem_dim(a).y],(void*)&(*ac)[yy + (xx+j)*lda], sizey*sizeof(T));
        }
    };
        
    template<typename T>
    struct reshape_l2r : public ambient::kernel_dispatch< reshape_l2r<T> > 
    {
        typedef void (*F)(const maquis::types::p_dense_matrix_impl<T>&, pinned maquis::types::p_dense_matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(const maquis::types::p_dense_matrix_impl<T>& left, pinned maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            ctxt_select("1 from ambient as reshape_l2r"); //if(!ctxt.involved()) return;
            block_2d_cycle_assign(ui_l_current(left)); 
            block_2d_cycle_pin(ui_l_current(right)); 
        }

        static inline void c(const maquis::types::p_dense_matrix_impl<T>& left, pinned maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        { // gs
            __A_TIME("ambient_reshape_l2r_c_kernel"); 
            __a_memptf<T, __a_memcpy>(right, dim2(0,0), right, dim2(0,0), dim2(ui_c_get_dim(right).x,ui_c_get_dim(right).y)); // refreshing updated memory
            for(size_t ss = 0; ss < sdim; ++ss){
                __a_memptf<T, __a_memcpy>(right, dim2(ss*rdim + right_offset, 0), 
                                          left,  dim2(0, ss*ldim + left_offset), 
                                          dim2( rdim, ldim ));
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct reshape_r2l : public ambient::kernel_dispatch< reshape_r2l<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
               const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& left, const maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            ctxt_select("1 from ambient as reshape_l2r"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(left)); 
            block_2d_cycle_assign(ui_l_current(right)); 
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& left, const maquis::types::p_dense_matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        { // gs
            __A_TIME("ambient_reshape_r2l_c_kernel"); 
            __a_memptf<T, __a_memcpy>(left, dim2(0,0), left, dim2(0,0), dim2(ui_c_get_dim(left).x,ui_c_get_dim(left).y)); // refreshing updated memory
            for(size_t ss = 0; ss < sdim; ++ss)
                __a_memptf<T, __a_memcpy>(left,  dim2(0, ss*ldim + left_offset), 
                                          right, dim2(ss*rdim + right_offset,0), 
                                          dim2( rdim, ldim ));
            __A_TIME_STOP
        }
    };
        
        
    template<typename T>
    struct rb_tensor_mpo : public ambient::kernel_dispatch< rb_tensor_mpo<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            ctxt_select("1 from ambient as rb_tensor_mpo"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(out)); 
            block_2d_cycle_assign(ui_l_current(in)); 
            block_2d_cycle_assign(ui_l_current(alfa)); 
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        { // gs
            __A_TIME("ambient_rb_tensor_mpo_c_kernel"); 
            __a_memptf<T, __a_memcpy>(out, dim2(0,0), out, dim2(0,0), dim2(ui_c_get_dim(out).x,ui_c_get_dim(out).y)); // refreshing updated memory
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T* alfad = ui_c_current(alfa)(ss2/ui_c_get_mem_dim(alfa).x, ss1/ui_c_get_mem_dim(alfa).y);
                    T  alfa_t = alfad[ss1%ui_c_get_mem_dim(alfa).y + ui_c_get_mem_dim(alfa).y*(ss2%ui_c_get_mem_dim(alfa).x)];
                    __a_memptf<T, __a_memscal>(out, dim2(out_offset + ss2*rdim, 0),
                                               in,  dim2(in_offset + ss1*rdim, 0),
                                               dim2(rdim, ldim), alfa_t);
                }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct lb_tensor_mpo : public ambient::kernel_dispatch< lb_tensor_mpo<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&,
               const size_t&, const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            ctxt_select("1 from ambient as rb_tensor_mpo"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(out)); 
            block_2d_cycle_assign(ui_l_current(in)); 
            block_2d_cycle_assign(ui_l_current(alfa)); 
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        { // gs
            __A_TIME("ambient_lb_tensor_mpo_c_kernel"); 
            __a_memptf<T, __a_memcpy>(out, dim2(0,0), out, dim2(0,0), dim2(ui_c_get_dim(out).x,ui_c_get_dim(out).y)); // refreshing updated memory
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T* alfad = ui_c_current(alfa)(ss2/ui_c_get_mem_dim(alfa).x, ss1/ui_c_get_mem_dim(alfa).y);
                    T  alfa_t = alfad[ss1%ui_c_get_mem_dim(alfa).y + ui_c_get_mem_dim(alfa).y*(ss2%ui_c_get_mem_dim(alfa).x)];
                    __a_memptf<T, __a_memscal>(out, dim2(0, out_offset + ss2*ldim),
                                               in,  dim2(0, in_offset + ss1*ldim),
                                               dim2(rdim, ldim), alfa_t);
                }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct scalar_norm : public ambient::kernel_dispatch< scalar_norm<T> > 
    {
        typedef void (*F)(pinned const maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, T*&);

        static inline void l(pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, T*& norm){
            ctxt_select("* from ambient as scalar_norm"); //if(!ctxt.involved()) return;
            block_outright_pin(ui_l_current(a));
        }

        static inline void c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, T*& norm){
            // gs
            __A_TIME("ambient_scalar_norm_c_kernel"); 
            T summ = 0;
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x,y);
            size_t lda = ui_c_get_mem_dim(a).y;
            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
            for(size_t i=0; i < sizey; i++)
                for(size_t j=0; j < sizex; j++)
                    summ += ad[i+j*lda]*ad[i+j*lda];
        
            *norm += summ;
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct scalar_overlap : public ambient::kernel_dispatch< scalar_overlap<T> > 
    {
        typedef void (*F)(pinned const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, T*&);

        static inline void l(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, const size_t& m, const size_t& n, T*& overlap){
            ctxt_select("* from ambient as scalar_overlap"); //if(!ctxt.involved()) return;
            block_outright_pin(ui_l_current(a));
            block_outright_assign(ui_l_current(b));
        }

        static inline void c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, const size_t& m, const size_t& n, T*& overlap){
            // gs
            __A_TIME("ambient_scalar_overlap_c_kernel"); 
            T summ = 0;
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x,y);
            T* bd = ui_c_current(b)(x,y);
            size_t lda = ui_c_get_mem_dim(a).y;
            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
            for(size_t i=0; i < sizey; i++)
                for(size_t j=0; j < sizex; j++)
                    summ += ad[i+j*lda]*bd[i+j*lda];
            *overlap += summ;
            __A_TIME_STOP
        }
    };
        
        
    template<typename T>
    struct add : public ambient::kernel_dispatch< add<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            ctxt_select("1 from ambient as add"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
            block_2d_cycle_assign(ui_l_current(b));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            // gs
            __A_TIME("ambient_add_c_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x,y);
            T* bd = ui_c_current(b)(x,y);
            T* ar = ui_c_updated(a)(x,y);
            size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
            for(size_t k = 0; k < size; k++)
                ar[k] = ad[k] + bd[k];
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct sub : public ambient::kernel_dispatch< sub<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            ctxt_select("1 from ambient as sub"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
            block_2d_cycle_assign(ui_l_current(b));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
            // gs
            __A_TIME("ambient_sub_c_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x,y);
            T* bd = ui_c_current(b)(x,y);
            T* ar = ui_c_updated(a)(x,y);
            size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
            for(size_t k = 0; k < size; k++)
                ar[k] = ad[k] + (-1)*bd[k];
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct scale : public ambient::kernel_dispatch< scale<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, const T*&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const T*& t){
            ctxt_select("1 from ambient as scale"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const T*& t){
            // gs
            __A_TIME("ambient_scale_c_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x, y);
            T* ar = ui_c_updated(a)(x, y);
        
            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
            size_t lda = ui_c_get_mem_dim(a).y;
            for(size_t j=0; j < sizex; j++)
            for(size_t i=0; i < sizey; i++)
            ar[j*lda+i] = ad[j*lda+i] * (*t);
            __A_TIME_STOP
        }
    };
        
    template<typename T, typename D>
    struct gemm_diagonal_lhs : public ambient::kernel_dispatch< gemm_diagonal_lhs<T,D> > 
    {
        typedef void (*F)(const maquis::types::p_dense_matrix_impl<D>&, pinned const maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&);

        static inline void l(const maquis::types::p_dense_matrix_impl<D>& a_diag, pinned const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c,
                      const size_t& m, const size_t& n, const size_t& k)
        {
            ctxt_select("1 from ambient as gemm_diagonal_lhs"); //if(!ctxt.involved()) return;
            block_2d_cycle_assign(ui_l_current(a_diag));
            block_2d_cycle_pin(ui_l_current(b));
            block_2d_cycle_assign(ui_l_current(c));
        }

        static inline void c(const maquis::types::p_dense_matrix_impl<D>& a_diag, pinned const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c,
                      const size_t& m, const size_t& n, const size_t& k){
            // gs
            __A_TIME("ambient_gemm_diagonal_lhs_c_kernel"); 
            size_t sizey = __a_get_limit_y(a_diag, m);
            int j = ctxt.get_block_id().y*ui_c_get_mem_dim(b).y;
            int size = ui_c_get_mem_dim(b).x;
            int lda  = sizeof(T)/sizeof(D)*ui_c_get_mem_dim(b).y;
            int ONE  = 1;
            D* bd = ui_c_current(b)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            D* cd = ui_c_updated(c)(ctxt.get_block_id().x, ctxt.get_block_id().y);
        
            for(int jj = 0 ; jj < sizey; jj++){
                 D* alpha = ui_c_current(a_diag)(0, (j+jj)/ui_c_get_mem_dim(a_diag).y);
        	     axpy(&size, &alpha[(j+jj)%ui_c_get_mem_dim(a_diag).y], &bd[jj], &lda, &cd[jj], &lda);
        	     if(sizeof(T) != sizeof(D)) axpy(&size, &alpha[(j+jj)%ui_c_get_mem_dim(a_diag).y], &bd[jj+1], &lda, &cd[jj], &lda); // for complex
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T, typename D>
    struct gemm_diagonal_rhs : public ambient::kernel_dispatch< gemm_diagonal_rhs<T,D> > 
    {
        typedef void (*F)(pinned const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<D>&, maquis::types::p_dense_matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&);

        static inline void l(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<D>& b_diag, maquis::types::p_dense_matrix_impl<T>& c,
                      const size_t& m, const size_t& n, const size_t& k)
        {
            ctxt_select("1 from ambient as gemm_diagonal_rhs"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
            block_2d_cycle_assign(ui_l_current(b_diag));
            block_2d_cycle_assign(ui_l_current(c));
        }

        static inline void c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<D>& b_diag, maquis::types::p_dense_matrix_impl<T>& c,
                      const size_t& m, const size_t& n, const size_t& k){
            // gs
            __A_TIME("ambient_gemm_diagonal_rhs_c_kernel"); 
            size_t sizex = std::min((ctxt.get_block_id().x+1)*ui_c_get_mem_dim(b_diag).y, n)-ctxt.get_block_id().x*ui_c_get_mem_dim(b_diag).y;
        
            int j = ctxt.get_block_id().x*ui_c_get_mem_dim(a).x;
            int size = sizeof(T)/sizeof(D)*ui_c_get_mem_dim(a).y; // for the case of complex
            int ONE = 1;
            D* ad = ui_c_current(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            D* cd = ui_c_updated(c)(ctxt.get_block_id().x, ctxt.get_block_id().y);
        
            for(int jj = 0 ; jj < sizex; jj++){
        	    D* alpha = ui_c_current(b_diag)(0, (j+jj)/ui_c_get_mem_dim(b_diag).y);
        	    axpy(&size, &alpha[(j+jj)%ui_c_get_mem_dim(b_diag).y], &ad[jj*ui_c_get_mem_dim(a).y], &ONE, &cd[jj*ui_c_get_mem_dim(c).y], &ONE);
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct trace : public ambient::kernel_dispatch< trace<T> > 
    {
        typedef void (*F)(pinned const maquis::types::p_dense_matrix_impl<T>&, const size_t&, T*&);

        static inline void l(pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& n, T*& trace){
            ctxt_select("* from ambient as trace"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a)); // we need only diagonal 
                                    // but we have to track the diagonal separately afterward
                                    // which is troublesome
        }

        static inline void c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& n, T*& trace){
            // gs
            __A_TIME("ambient_trace_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            size_t ld = ui_c_get_mem_dim(a).y;
            size_t sd = ui_c_get_mem_dim(a).x;
            T* ad = ui_c_current(a)(x, y);
        
            if((y+1)*ld <= x*sd) return;
            if(y*ld >= (x+1)*sd) return;
            size_t sizex = std::min(n,(x+1)*sd);
            for(size_t jj = x*sd; jj < sizex; jj++){
                if(y*ld > jj) continue;
                if((y+1)*ld <= jj) continue;
               *trace += ad[jj % ld + (jj%sd)*ld];
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct transpose_out : public ambient::kernel_dispatch< transpose_out<T> > 
    { // only for square blocks
        typedef void (*F)(pinned const maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&);

        static inline void l(pinned const maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<T>& t, const size_t& m, const size_t& n){
            ctxt_select("1 from ambient as transpose_out_l"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
            block_2d_cycle_assign(ui_l_current(t));
        }

        static inline void c(pinned const maquis::types::p_dense_matrix_impl<T>& a, maquis::types::p_dense_matrix_impl<T>& t, const size_t& m, const size_t& n){
            // gs
            __A_TIME("ambient_transpose_out_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            T* od = ui_c_current(a)(x,y);

            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
            size_t mlda = ui_c_get_mem_dim(a).y;
            size_t tlda = ui_c_get_mem_dim(t).y;

            size_t txx = y*ui_c_get_mem_dim(a).y;
            size_t tyy = x*ui_c_get_mem_dim(a).x;

            for(size_t j=0; j < sizex; ++j){
                for(size_t i = 0; i < sizey; ++i){
                    T* td = ui_c_updated(t)((txx+i) / ui_c_get_mem_dim(t).x, (tyy+j) / ui_c_get_mem_dim(t).y);
                    size_t ti = (tyy+j) % ui_c_get_mem_dim(t).y;
                    size_t tj = (txx+i) % ui_c_get_mem_dim(t).x;
                    td[ti+tj*tlda] = od[i+j*mlda];
                }
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct init_value : public ambient::kernel_dispatch< init_value<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&, const T&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const T& value){
            ctxt_select("1 from ambient as init_value"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const T& value){
            __A_TIME("ambient_init_value_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
        
            T* ad = ui_c_updated(a)(x,y);
            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
        
            for(size_t j=0; j < sizex; ++j){
                for(size_t i = 0; i < sizey; ++i){
                    ad[i+j*ui_c_get_mem_dim(a).y] = value; // not a memset due to complex
                }
            }
            __A_TIME_STOP
        }
    };
       
    template<typename T>
    struct init_random : public ambient::kernel_dispatch< init_random<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&);
     
        template<typename T> static inline void randomize(T* ad){ *ad = drand48(); }
        template<typename T> static inline void randomize(std::complex<T>* ad){
            (*ad).real() = drand48();
            (*ad).imag() = drand48();
        }

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
            ctxt_select("1 from ambient as init_random"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }
        
        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
            __A_TIME("ambient_init_random_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            size_t ld = ui_c_get_mem_dim(a).y;
            size_t sd = ui_c_get_mem_dim(a).x;
            T* ad = ui_c_updated(a)(x,y);
            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
          
            for(size_t jj = 0; jj < sizex; jj++){
                for(size_t ii = 0; ii < sizey; ii++)
                    randomize((ad+(jj*ld+ii)));
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct init_identity : public ambient::kernel_dispatch< init_identity<T> > 
    {
        typedef void (*F)(pinned maquis::types::p_dense_matrix_impl<T>&, const size_t&, const size_t&);

        static inline void l(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
            ctxt_select("1 from ambient as init_identity"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a));
        }

        static inline void c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
            __A_TIME("ambient_init_identity_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            size_t ld = ui_c_get_mem_dim(a).y;
            size_t sd = ui_c_get_mem_dim(a).x;
            T* ad = ui_c_updated(a)(x,y);
            if((y+1)*ld <= x*sd) return;
            if(y*ld >= (x+1)*sd) return;
            size_t sizex = std::min(m,n); // respecting borders
            sizex = std::min(sizex,(x+1)*sd);
            for(size_t jj = x*sd; jj < sizex; jj++){
                if(y*ld > jj) continue;
                if((y+1)*ld <= jj) continue;
                ad[jj % ld + (jj%sd)*ld] = 1.;
            }
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct validation : public ambient::kernel_dispatch< validation<T> > 
    {
        typedef void (*F)(pinned const maquis::types::p_dense_matrix_impl<T>&, const maquis::types::p_dense_matrix_impl<T>&, int*&);

        static inline void l(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, int*& ret){
            ctxt_select("1 from ambient as validation"); //if(!ctxt.involved()) return;
            block_2d_cycle_pin(ui_l_current(a)); 
            block_2d_cycle_assign(ui_l_current(b)); 
        }
        
        static inline void c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, int*& ret){ // see paper for Reference Dongara 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x, y); 
            T* bd = ui_c_current(b)(x, y); 
            double res(0.0); 
            double epsilon = std::numeric_limits<double>::epsilon();
            size_t position_x(0),position_y(0),position_xy(0); 
       
            for(size_t ii=0; ii < ui_c_get_mem_dim(a).y; ++ii){
                for(size_t jj=0; jj < ui_c_get_mem_dim(a).x; ++jj){
                    position_x = x*ui_c_get_mem_dim(a).x+jj;
                    position_y = y*ui_c_get_mem_dim(a).y+ii;
                    if(position_x < std::min(ui_c_get_dim(a).x,ui_c_get_dim(b).x) && 
                       position_y < std::min(ui_c_get_dim(a).y,ui_c_get_dim(b).y)  )
                    {
                        position_xy = jj*ui_c_get_mem_dim(a).y+ii;
                        res = (norm(ad[position_xy])-norm(bd[position_xy]))/fabs(epsilon*norm(bd[position_xy])); // to do : rotation pb  with complex to change
                        if(res > 256){ // 16 is recommended by Dongara, 256 because lapack gives != runs after runs
                            std::cout << position_y << " " << position_x << " : " << ad[position_xy] << " " << bd[position_xy] << std::endl;
                            *ret = 0; // test failed return 0 (bool false)
                        }
                    }
                }
            }
        }
    };

    // {{{ MKL LAPACK kernels

    template<typename T>
    struct svd : public ambient::kernel_dispatch< svd<T> > 
    {
        typedef void (*F)(const maquis::types::p_dense_matrix_impl<T>&, int&, int&, maquis::types::p_dense_matrix_impl<T>&, 
                          maquis::types::p_dense_matrix_impl<T>&, maquis::types::p_dense_matrix_impl<double>&);

        static inline void l(const maquis::types::p_dense_matrix_impl<T>& a, int& m, int& n, maquis::types::p_dense_matrix_impl<T>& u, 
                      maquis::types::p_dense_matrix_impl<T>& vt, maquis::types::p_dense_matrix_impl<double>& s)
        {
            ctxt_select("1 from ambient as svd"); //if(!ctxt.involved()) return;
            block_outright_conditional_assign(ui_l_current(s));
            block_2d_cycle_conditional_assign(ui_l_current(a));
            block_2d_cycle_conditional_assign(ui_l_current(u));
            block_2d_cycle_conditional_assign(ui_l_current(vt));
        }

        static inline void c(const maquis::types::p_dense_matrix_impl<T>& a, int& m, int& n, maquis::types::p_dense_matrix_impl<T>& u, 
                     maquis::types::p_dense_matrix_impl<T>& vt, maquis::types::p_dense_matrix_impl<double>& s)
        {
            // gs
            __A_TIME("ambient_svd_c_kernel"); 
        /* Locals */
            int lda = ui_c_get_grid_dim(a).y*ui_c_get_mem_dim(a).y;
            int ldu = ui_c_get_grid_dim(u).y*ui_c_get_mem_dim(u).y;
            int ldvt = ui_c_get_grid_dim(vt).y*ui_c_get_mem_dim(vt).y;
            int info, lwork;
            T wkopt;
            T* work;
            double* rwork = new double[5*std::min(lda,ldu)]; // C - useless for double but need for complex 
            T* ad;
            T* ud;
            T* vtd;
            double* sd;
        
            if(ui_c_get_grid_dim(a) == 1){
                ad  = ui_c_current(a) (0,0);
                ud  = ui_c_updated(u) (0,0);
                vtd = ui_c_updated(vt)(0,0);
                sd  = ui_c_updated(s) (0,0);
            }else{
                ad  = (T*)__a_solidify(a);
                ud  = (T*)__a_solidify(u);
                vtd = (T*)__a_solidify(vt);
                sd  = (double*)__a_solidify(s);
            }
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
            if(ui_c_get_grid_dim(a) != 1){
                __a_disperse(ud, u);
                __a_disperse(vtd, vt);
                __a_disperse(sd, s);
            }
            free(work);
            __A_TIME_STOP
        }
    };
        
    template<typename T>
    struct syev : public ambient::kernel_dispatch< syev<T> > 
    {
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, int&, maquis::types::p_dense_matrix_impl<T>&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& a, int& m, maquis::types::p_dense_matrix_impl<T>& w){
            ctxt_select("1 from ambient as syev"); //if(!ctxt.involved()) return;
            block_2d_cycle_conditional_assign(ui_l_current(a));
            block_2d_cycle_conditional_assign(ui_l_current(w));
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& a, int& m, maquis::types::p_dense_matrix_impl<T>& w){
            int lda = ui_c_get_grid_dim(a).y*ui_c_get_mem_dim(a).y;
            int info, lwork = -1;
            double wkopt;
            double* work;
            double* ad = (double*)__a_solidify(a);
            double* wd = ui_c_get_grid_dim(a) == 1 ? ui_c_updated(w)(0,0) : (double*)__a_solidify(w);
             
            dsyev_("V","U",&m,ad,&lda,wd,&wkopt,&lwork,&info);
            lwork = (int)wkopt;
            work = (double*)malloc( lwork*sizeof(double) );
            dsyev_("V","U",&m,ad,&lda,wd,work,&lwork,&info);
        
            if( info > 0 ) {
                printf( "The algorithm computing SYEV failed to converge.\n" );
                exit( 1 );
            }
        
            __a_disperse(ad, a);
            if(ui_c_get_grid_dim(a) != 1) 
                __a_disperse(wd, w);
            free(work); 
        }
    };
        
    template<typename T>
    struct heev : public ambient::kernel_dispatch< heev<T> > 
    {
        typedef void (*F)(maquis::types::p_dense_matrix_impl<T>&, const size_t&, maquis::types::p_dense_matrix_impl<double>&);

        static inline void l(maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, maquis::types::p_dense_matrix_impl<double>& w){
            ctxt_select("1 from ambient as heev"); //if(!ctxt.involved()) return;
            block_2d_cycle_conditional_assign(ui_l_current(a));
            block_2d_cycle_conditional_assign(ui_l_current(w)); // C - block_outright(w) is possible, if yes remove solidify and disperse for w 
        }

        static inline void c(maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, maquis::types::p_dense_matrix_impl<double>& w){
            // gs
            __A_TIME("ambient_heev_c_kernel"); 
            int lda = ui_c_get_grid_dim(a).y*ui_c_get_mem_dim(a).y;
            int info, lwork = -1;
        
            double wkopt;
            double* work;
            int am = (int)m; // for mkl (int*)
        
            double* ad = (double*)__a_solidify(a);
            double* wd = (double*)__a_solidify(w);
       
            if(ui_c_current(a).get_layout().mem_dim != ui_c_updated(a).get_layout().mem_dim) printf("Dimensions don't match!\n"); 

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
         
            __a_disperse(ad, a);
            __a_disperse(wd, w);
            free(work);
            __A_TIME_STOP
        }
    };

    // }}}

}
#endif
