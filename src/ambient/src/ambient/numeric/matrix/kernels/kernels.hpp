#ifndef __AMBIENT_NUMERIC_MATRIX_KERNELS_HPP__
#define __AMBIENT_NUMERIC_MATRIX_KERNELS_HPP__

namespace ambient { namespace numeric { namespace kernels {

    using ambient::numeric::matrix_impl;

    template<typename T>
    struct gemm_inplace : public kernel< gemm_inplace<T> > 
    {
        typedef void(gemm_inplace::*F)(matrix_impl<T>&, 
                         const matrix_impl<T>&);

        inline void l(matrix_impl<T>& a, const matrix_impl<T>& b){
            this->ctxt_select("1 from ambient as gemm_inplace"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
            this->block_2d_cycle_assign(b);
        }

        inline void c(matrix_impl<T>& a, const matrix_impl<T>& b){
            __A_TIME_C("ambient_gemm_inplace_c_kernel");
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
                        //if(!ui_p_updated(a)(y,*zy).trylock()){ zy++; continue; }
                        T* ad = ui_c_current(a)(x,*zy);
                        T* cd = ui_p_updated(a)(y,*zy); // a(x,z) x b(y,x) => c(y,z)
                        gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
                        //ui_p_updated(a)(y,*zy).unlock();
                        L.erase(zy++);
                    }
                }
                y += ui_c_get_grid_dim(a).y;
            }
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct gemm_general : public kernel< gemm_general<T> > 
    {
        typedef void(gemm_general::*F)(const matrix_impl<T>&, 
                         const matrix_impl<T>&, 
                               matrix_impl<T>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<T>& b, matrix_impl<T>& c){
            this->ctxt_select("1 from ambient as gemm"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
            this->block_2d_cycle_assign(b);
            this->block_2d_cycle_assign(c);
        }

        inline void c(const matrix_impl<T>& a, const matrix_impl<T>& b, matrix_impl<T>& c){
            // gs
            __A_TIME_C("ambient_gemm_c_kernel");
            //if(ui_c_get_grid_dim(a) == 1 && ui_c_get_grid_dim(b) == 1) // use gemm atomic for this case
            if(ui_c_get_mem_dim(a) != ui_c_get_mem_dim(b) || ui_c_get_mem_dim(a) != ui_c_get_mem_dim(c)){
                T* ad = (T*)__a_solidify<T>(a);
                T* bd = (T*)__a_solidify<T>(b);
                T* cd = (T*)__a_solidify<T>(c);
                int m   = ui_c_get_dim(a).y;
                int n   = ui_c_get_dim(b).x;
                int k   = ui_c_get_dim(b).y;
                int lda = ui_c_get_mem_dim(a).y*ui_c_get_grid_dim(a).y;
                int ldb = ui_c_get_mem_dim(b).y*ui_c_get_grid_dim(b).y;
                int ldc = ui_c_get_mem_dim(c).y*ui_c_get_grid_dim(c).y;
                T alpha(1.0); 
                T beta(0.0);
                gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
                __a_disperse<T>(cd, c);
                __A_TIME_C_STOP
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
                        //if(!ui_p_updated(c)(y,*zy).trylock()){ zy++; continue; }
                        T* ad = ui_c_current(a)(x,*zy);
                        T* cd = ui_p_updated(c)(y,*zy); // a(x,z) x b(y,x) => c(y,z)
                        gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
                        //ui_p_updated(c)(y,*zy).unlock();
                        L.erase(zy++);
                    }
                }
                y += a_grid_dim.y;
            }
            __A_TIME_C_STOP
        }
    };


    template<typename T>
    struct copy : public kernel< copy<T> > 
    {
        typedef void(copy::*F)(matrix_impl<T>&, const matrix_impl<T>&);

        inline void l(matrix_impl<T>& ac, const matrix_impl<T>& a){
            this->ctxt_select("1 from ambient as copy"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
            this->block_2d_cycle_assign(ac);
        }

        inline void c(matrix_impl<T>& ac, const matrix_impl<T>& a){
            // gs
            __A_TIME_C("ambient_copy_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            T* a_elements  = ui_c_current(a)(x,y);
            T* ac_elements = ui_w_updated(ac)(x,y);
            memcpy(ac_elements, a_elements, sizeof(T)*ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y);
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct remove_rows : public kernel< remove_rows<T> > 
    {
        typedef void (remove_rows::*F)(matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& a, const size_t& i_mark, const size_t& k){
            this->ctxt_select("1 from ambient as remove_rows"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }
        
        inline void c(matrix_impl<T>& a, const size_t& i_mark, const size_t& k){
            __A_TIME_C("ambient_remove_rows_c_kernel"); 
            size_t numcols = ui_c_get_dim(a).x;
            size_t numrows = ui_c_get_dim(a).y;
            if((T*)ui_c_current(a)(0,0) != (T*)ui_w_updated(a)(0,0)) 
                __a_memptf_reverse<T, __a_memcpy>(a, dim2(0,0), a, dim2(0,0), dim2(numcols, i_mark));
            __a_memptf_reverse<T, __a_memcpy>(a, dim2(0,i_mark), a, dim2(0,k+i_mark), dim2(numcols,numrows-k-i_mark));
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct remove_cols : public kernel< remove_cols<T> > 
    {
        typedef void (remove_cols::*F)(matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& a, const size_t& j_mark, const size_t& k){
            this->ctxt_select("1 from ambient as remove_cols"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }

        inline void c(matrix_impl<T>& a, const size_t& j_mark, const size_t& k){
            __A_TIME_C("ambient_remove_cols_c_kernel"); 
            size_t numcols = ui_c_get_dim(a).x;
            size_t numrows = ui_c_get_dim(a).y;
            if((T*)ui_c_current(a)(0,0) != (T*)ui_w_updated(a)(0,0)) 
                __a_memptf_reverse<T, __a_memcpy>(a, dim2(0,0), a, dim2(0,0), dim2(j_mark, numrows));
            __a_memptf_reverse<T, __a_memcpy>(a, dim2(j_mark,0), a, dim2(k+j_mark,0), dim2(numcols-k-j_mark,numrows));
            __A_TIME_C_STOP
        }
    };

    template<typename T>
    struct resize : public kernel< resize<T> > 
    {
        typedef void (resize::*F)(matrix_impl<T>&, const matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& r, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            this->ctxt_select("1 from ambient as resize"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
            this->block_2d_cycle_assign(r);
        }

        inline void c(matrix_impl<T>& r, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            __A_TIME_C("ambient_resize_c_kernel"); 
            __a_memptf_reverse<T, __a_memcpy>(r, dim2(0,0), a, dim2(0,0), dim2(n,m));
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct sqrt_diagonal : public kernel< sqrt_diagonal<T> > 
    {
        typedef void (sqrt_diagonal::*F)(matrix_impl<T>&);

        inline void l(matrix_impl<T>& a){
            this->ctxt_select("1 from ambient as sqrt_diagonal"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }

        inline void c(matrix_impl<T>& a){
            __A_TIME_C("ambient_sqrt_diagonal_c_kernel"); 
            T* ad = ui_c_current(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            T* sd = ui_w_updated(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            size_t size = ui_c_get_mem_dim(a).y;
            for(int i=0; i < size; i++)
                sd[i] = std::sqrt(ad[i]);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct exp_diagonal : public kernel< exp_diagonal<T> > 
    {
        typedef void (exp_diagonal::*F)(matrix_impl<T>&, const T&);

        inline void l(matrix_impl<T>& a, const T& alfa){
            this->ctxt_select("1 from ambient as exp_diagonal"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }

        inline void c(matrix_impl<T>& a, const T& alfa){
            __A_TIME_C("ambient_exp_diagonal_c_kernel"); 
            T* ad = ui_c_current(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            T* sd = ui_w_updated(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            size_t size = ui_c_get_mem_dim(a).y;
            for(int i=0; i < size; i++)
                sd[i] = std::exp(ad[i]*alfa);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct exp_diagonal_rc : public kernel< exp_diagonal_rc<T> > 
    {
        typedef void (exp_diagonal_rc::*F)(matrix_impl< std::complex<T> >&, const matrix_impl<T>&, const std::complex<T>&);

        inline void l(matrix_impl< std::complex<T> >& e, const matrix_impl<T>& a, const std::complex<T>& alfa){
            this->ctxt_select("1 from ambient as exp_diagonal"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(e);
            this->block_2d_cycle_assign(a);
        }

        inline void c(matrix_impl< std::complex<T> >& e, const matrix_impl<T>& a, const std::complex<T>& alfa){
            __A_TIME_C("ambient_exp_diagonal_rc_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            std::complex<T>* ed = ui_w_updated(e)(x,y);
            T* ad = ui_c_current(a)(x,y);
            size_t size = ui_c_get_mem_dim(e).y;
            for(int i=0; i < size; i++)
                ed[i] = std::exp(ad[i]*alfa);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct push_back_sqr_gt : public kernel< push_back_sqr_gt<T> > 
    {
        typedef void (push_back_sqr_gt::*F)(const matrix_impl<T>&, std::vector<T>*&);

        inline void l(const matrix_impl<T>& a, std::vector<T>*& ac){
            this->ctxt_select("* from ambient as push_back_sqr_gt"); //if(!ctxt.involved()) return;
            this->block_outright_pin(a);
        }

        inline void c(const matrix_impl<T>& a, std::vector<T>*& ac){
            // gs
            __A_TIME_C("ambient_push_back_sqr_gt_c_kernel"); 
            T* ad = ui_c_current(a)(ctxt.get_block_id().x, ctxt.get_block_id().y);
            size_t sizey = __a_get_limit_y(a);
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
        typedef void (cast_to_vector::*F)(std::vector<T>*&, const matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(std::vector<T>*& ac, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            this->ctxt_select("* from ambient as cast_to_vector"); //if(!ctxt.involved()) return;
            this->block_outright_pin(a);
        }

        inline void c(std::vector<T>*& ac, const matrix_impl<T>& a, const size_t& m, const size_t& n){
            // gs
            __A_TIME_C("ambient_cast_to_dense_c_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            int xx = ui_c_get_mem_dim(a).x*x; 
            int yy = ui_c_get_mem_dim(a).y*y; // conversion cartersia coordinates alps / ambient
            size_t offset;
            T* ad = ui_c_current(a)(x,y);
           
            size_t sizex = __a_get_limit_x(a, n);
            size_t sizey = __a_get_limit_y(a, m);
            size_t lda = ui_c_get_mem_dim(a).y;
        
            for(int j=0; j < sizex; ++j){
                offset = yy + (xx+j)*m;
                memcpy((void*)&(*ac)[offset],(void*)&ad[j*lda], sizey*sizeof(T));  
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct cast_from_vector : public kernel< cast_from_vector<T> > 
    {
        typedef void (cast_from_vector::*F)(const std::vector<T>*&, matrix_impl<T>&, const size_t&, const size_t&, const size_t&);

        inline void l(const std::vector<T>*& ac, matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& lda){
            this->ctxt_select("1 from ambient as cast_from_vector"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }

        inline void c(const std::vector<T>*& ac, matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& lda){
            __A_TIME_C("ambient_cast_from_vector_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            size_t xx = ui_c_get_mem_dim(a).x*x;
            size_t yy = ui_c_get_mem_dim(a).y*y;
            size_t sizex = __a_get_limit_x(a, n);
            size_t sizey = __a_get_limit_y(a, m);

            T* ad = ui_w_updated(a)(x,y);
            for(int j=0; j < sizex; ++j)
                memcpy((void*)&ad[j*ui_c_get_mem_dim(a).y],(void*)&(*ac)[yy + (xx+j)*lda], sizey*sizeof(T));
            __A_TIME_C_STOP 
        }
    };
        
    template<typename T>
    struct reshape_l2r : public kernel< reshape_l2r<T> > 
    {
        typedef void (reshape_l2r::*F)(const matrix_impl<T>&, matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(const matrix_impl<T>& left, matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as reshape_l2r"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_assign(left); 
            this->block_2d_cycle_pin(right); 
        }

        inline void c(const matrix_impl<T>& left, matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        { // gs
            __A_TIME_C("ambient_reshape_l2r_c_kernel"); 
            if((T*)ui_c_current(right)(0,0) != (T*)ui_w_updated(right)(0,0)) 
                __a_memptf<T, __a_memcpy>(right, dim2(0,0), right, dim2(0,0), ui_c_get_dim(right)); // refreshing updated memory
            for(size_t ss = 0; ss < sdim; ++ss){
                __a_memptf<T, __a_memcpy>(right, dim2(ss*rdim + right_offset, 0), 
                                          left,  dim2(0, ss*ldim + left_offset), 
                                          dim2( rdim, ldim ));
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct reshape_r2l : public kernel< reshape_r2l<T> > 
    {
        typedef void (reshape_r2l::*F)(matrix_impl<T>&, const matrix_impl<T>&,
               const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& left, const matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as reshape_l2r"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(left); 
            this->block_2d_cycle_assign(right); 
        }

        inline void c(matrix_impl<T>& left, const matrix_impl<T>& right,
                      const size_t& left_offset, const size_t& right_offset, 
                      const size_t& sdim, const size_t& ldim, const size_t& rdim)
        { // gs
            __A_TIME_C("ambient_reshape_r2l_c_kernel"); 
            if((T*)ui_c_current(left)(0,0) != (T*)ui_w_updated(left)(0,0)) 
                __a_memptf<T, __a_memcpy>(left, dim2(0,0), left, dim2(0,0), ui_c_get_dim(left)); // refreshing updated memory
            for(size_t ss = 0; ss < sdim; ++ss)
                __a_memptf<T, __a_memcpy>(left,  dim2(0, ss*ldim + left_offset), 
                                          right, dim2(ss*rdim + right_offset,0), 
                                          dim2( rdim, ldim ));
            __A_TIME_C_STOP
        }
    };
        
        
    template<typename T>
    struct rb_tensor_mpo : public kernel< rb_tensor_mpo<T> > 
    {
        typedef void (rb_tensor_mpo::*F)(matrix_impl<T>&, const matrix_impl<T>&, const matrix_impl<T>&,
                          const size_t&, const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as rb_tensor_mpo"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(out); 
            this->block_2d_cycle_assign(in); 
            this->block_2d_cycle_assign(alfa); 
        }

        inline void c(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        { // gs
            __A_TIME_C("ambient_rb_tensor_mpo_c_kernel"); 
            if((T*)ui_c_current(out)(0,0) != (T*)ui_w_updated(out)(0,0)) 
                __a_memptf<T, __a_memcpy>(out, dim2(0,0), out, dim2(0,0), ui_c_get_dim(out)); // refreshing updated memory
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T* alfad = ui_c_current(alfa)(ss2/ui_c_get_mem_dim(alfa).x, ss1/ui_c_get_mem_dim(alfa).y);
                    T  alfa_t = alfad[ss1%ui_c_get_mem_dim(alfa).y + ui_c_get_mem_dim(alfa).y*(ss2%ui_c_get_mem_dim(alfa).x)];
                    __a_memptf<T, __a_memscal>(out, dim2(out_offset + ss2*rdim, 0),
                                               in,  dim2(in_offset + ss1*rdim, 0),
                                               dim2(rdim, ldim), alfa_t);
                }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct lb_tensor_mpo : public kernel< lb_tensor_mpo<T> > 
    {
        typedef void (lb_tensor_mpo::*F)(matrix_impl<T>&, const matrix_impl<T>&, const matrix_impl<T>&,
               const size_t&, const size_t&, const size_t&, const size_t&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        {
            this->ctxt_select("1 from ambient as rb_tensor_mpo"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(out); 
            this->block_2d_cycle_assign(in); 
            this->block_2d_cycle_assign(alfa); 
        }

        inline void c(matrix_impl<T>& out, const matrix_impl<T>& in, const matrix_impl<T>& alfa,
                      const size_t& out_offset, const size_t& in_offset, 
                      const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
        { // gs
            __A_TIME_C("ambient_lb_tensor_mpo_c_kernel"); 
            if((T*)ui_c_current(out)(0,0) != (T*)ui_w_updated(out)(0,0)) 
                __a_memptf<T, __a_memcpy>(out, dim2(0,0), out, dim2(0,0), ui_c_get_dim(out)); // refreshing updated memory
            for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
                for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                    T* alfad = ui_c_current(alfa)(ss2/ui_c_get_mem_dim(alfa).x, ss1/ui_c_get_mem_dim(alfa).y);
                    T  alfa_t = alfad[ss1%ui_c_get_mem_dim(alfa).y + ui_c_get_mem_dim(alfa).y*(ss2%ui_c_get_mem_dim(alfa).x)];
                    __a_memptf<T, __a_memscal>(out, dim2(0, out_offset + ss2*ldim),
                                               in,  dim2(0, in_offset + ss1*ldim),
                                               dim2(rdim, ldim), alfa_t);
                }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scalar_norm : public kernel< scalar_norm<T> > 
    {
        typedef void (scalar_norm::*F)(const matrix_impl<T>&, const size_t&, const size_t&, future<double>&);

        inline void l(const matrix_impl<T>& a, const size_t& m, const size_t& n, future<double>& norm){
            this->ctxt_select("* from ambient as scalar_norm"); //if(!ctxt.involved()) return;
            this->block_outright_pin(a);
        }

        inline void c(const matrix_impl<T>& a, const size_t& m, const size_t& n, future<double>& norm){
            // gs
            __A_TIME_C("ambient_scalar_norm_c_kernel"); 
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
        
            norm.get_value() += alps::numeric::real(summ);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct overlap : public kernel< overlap<T> > 
    {
        typedef void (overlap::*F)(const matrix_impl<T>&, const matrix_impl<T>&, const size_t&, const size_t&, future<T>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<T>& b, const size_t& m, const size_t& n, future<T>& overlap){
            this->ctxt_select("* from ambient as overlap"); //if(!ctxt.involved()) return;
            this->block_outright_pin(a);
            this->block_outright_assign(b);
        }

        inline void c(const matrix_impl<T>& a, const matrix_impl<T>& b, const size_t& m, const size_t& n, future<T>& overlap){
            // gs
            __A_TIME_C("ambient_scalar_overlap_c_kernel"); 
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
            overlap.get_value() += summ;
            __A_TIME_C_STOP
        }
    };
        
        
    template<typename T>
    struct add : public kernel< add<T> > 
    {
        typedef void (add::*F)(matrix_impl<T>&, const matrix_impl<T>&);

        inline void l(matrix_impl<T>& a, const matrix_impl<T>& b){
            this->ctxt_select("1 from ambient as add"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
            this->block_2d_cycle_assign(b);
        }

        inline void c(matrix_impl<T>& a, const matrix_impl<T>& b){
            // gs
            __A_TIME_C("ambient_add_c_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x,y);
            T* bd = ui_c_current(b)(x,y);
            T* ar = ui_w_updated(a)(x,y);
            size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
            for(size_t k = 0; k < size; k++)
                ar[k] = ad[k] + bd[k];
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct sub : public kernel< sub<T> > 
    {
        typedef void (sub::*F)(matrix_impl<T>&, const matrix_impl<T>&);

        inline void l(matrix_impl<T>& a, const matrix_impl<T>& b){
            this->ctxt_select("1 from ambient as sub"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
            this->block_2d_cycle_assign(b);
        }

        inline void c(matrix_impl<T>& a, const matrix_impl<T>& b){
            // gs
            __A_TIME_C("ambient_sub_c_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x,y);
            T* bd = ui_c_current(b)(x,y);
            T* ar = ui_w_updated(a)(x,y);
            size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
            for(size_t k = 0; k < size; k++)
                ar[k] = ad[k] + (-1)*bd[k];
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct scale : public kernel< scale<T> > 
    {
        typedef void (scale::*F)(matrix_impl<T>&, const size_t&, const size_t&, const future<T>&);

        inline void l(matrix_impl<T>& a, const size_t& m, const size_t& n, const future<T>& t){
            this->ctxt_select("1 from ambient as scale"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }

        inline void c(matrix_impl<T>& a, const size_t& m, const size_t& n, const future<T>& t){
            // gs
            __A_TIME_C("ambient_scale_c_kernel"); 
            int x = ctxt.get_block_id().x;
            int y = ctxt.get_block_id().y;
            T* ad = ui_c_current(a)(x, y);
            T* ar = ui_w_updated(a)(x, y);
        
            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
            size_t lda = ui_c_get_mem_dim(a).y;
            for(size_t j=0; j < sizex; j++)
            for(size_t i=0; i < sizey; i++)
            ar[j*lda+i] = ad[j*lda+i] * t.get_value();
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct trace : public kernel< trace<T> > 
    {
        typedef void (trace::*F)(const matrix_impl<T>&, const size_t&, future<T>&);

        inline void l(const matrix_impl<T>& a, const size_t& n, future<T>& trace){
            this->ctxt_select("* from ambient as trace"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);    // we need only diagonal 
                                            // but we have to track the diagonal separately afterward
                                            // which is troublesome
        }

        inline void c(const matrix_impl<T>& a, const size_t& n, future<T>& trace){
            // gs
            __A_TIME_C("ambient_trace_c_kernel"); 
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
                trace.get_value() += ad[jj % ld + (jj%sd)*ld];
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct transpose_out : public kernel< transpose_out<T> > 
    { // only for square blocks
        typedef void (transpose_out::*F)(const matrix_impl<T>&, matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(const matrix_impl<T>& a, matrix_impl<T>& t, const size_t& m, const size_t& n){
            this->ctxt_select("1 from ambient as transpose_out_l"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
            this->block_2d_cycle_assign(t);
        }

        inline void c(const matrix_impl<T>& a, matrix_impl<T>& t, const size_t& m, const size_t& n){
            // gs
            __A_TIME_C("ambient_transpose_out_c_kernel"); 
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
                    T* td = ui_w_updated(t)((txx+i) / ui_c_get_mem_dim(t).x, (tyy+j) / ui_c_get_mem_dim(t).y);
                    size_t ti = (tyy+j) % ui_c_get_mem_dim(t).y;
                    size_t tj = (txx+i) % ui_c_get_mem_dim(t).x;
                    td[ti+tj*tlda] = od[i+j*mlda];
                }
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct init_value : public kernel< init_value<T> > 
    {
        typedef void (init_value::*F)(matrix_impl<T>&, const size_t&, const size_t&, const T&);

        inline void l(matrix_impl<T>& a, const size_t& m, const size_t& n, const T& value){
            this->ctxt_select("1 from ambient as init_value"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }

        inline void c(matrix_impl<T>& a, const size_t& m, const size_t& n, const T& value){
            __A_TIME_C("ambient_init_value_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
        
            T* ad = ui_w_updated(a)(x,y);
            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
        
            for(size_t j=0; j < sizex; ++j){
                for(size_t i = 0; i < sizey; ++i){
                    ad[i+j*ui_c_get_mem_dim(a).y] = value; // not a memset due to complex
                }
            }
            __A_TIME_C_STOP
        }
    };
       
    template<typename T>
    struct init_random : public kernel< init_random<T> > 
    {
        typedef void (init_random::*F)(matrix_impl<T>&, const size_t&, const size_t&);
     
        template<typename T> inline void randomize(T* ad){ *ad = drand48(); }
        template<typename T> inline void randomize(std::complex<T>* ad){
            ad->real(drand48());
            ad->imag(drand48());
        }

        inline void l(matrix_impl<T>& a, const size_t& m, const size_t& n){
            this->ctxt_select("1 from ambient as init_random"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }
        
        inline void c(matrix_impl<T>& a, const size_t& m, const size_t& n){
            __A_TIME_C("ambient_init_random_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            size_t ld = ui_c_get_mem_dim(a).y;
            size_t sd = ui_c_get_mem_dim(a).x;
            T* ad = ui_w_updated(a)(x,y);
            size_t sizey = __a_get_limit_y(a, m);
            size_t sizex = __a_get_limit_x(a, n);
          
            for(size_t jj = 0; jj < sizex; jj++){
                for(size_t ii = 0; ii < sizey; ii++)
                    randomize((ad+(jj*ld+ii)));
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct init_identity : public kernel< init_identity<T> > 
    {
        typedef void (init_identity::*F)(matrix_impl<T>&, const size_t&, const size_t&);

        inline void l(matrix_impl<T>& a, const size_t& m, const size_t& n){
            this->ctxt_select("1 from ambient as init_identity"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a);
        }

        inline void c(matrix_impl<T>& a, const size_t& m, const size_t& n){
            __A_TIME_C("ambient_init_identity_c_kernel"); 
            size_t x = ctxt.get_block_id().x;
            size_t y = ctxt.get_block_id().y;
            size_t ld = ui_c_get_mem_dim(a).y;
            size_t sd = ui_c_get_mem_dim(a).x;
            T* ad = ui_r_updated(a)(x,y);
            if((y+1)*ld <= x*sd) return;
            if(y*ld >= (x+1)*sd) return;
            size_t sizex = std::min(m,n); // respecting borders
            sizex = std::min(sizex,(x+1)*sd);
            for(size_t jj = x*sd; jj < sizex; jj++){
                if(y*ld > jj) continue;
                if((y+1)*ld <= jj) continue;
                ad[jj % ld + (jj%sd)*ld] = 1.;
            }
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct validation : public kernel< validation<T> > 
    {
        typedef void (validation::*F)(const matrix_impl<T>&, const matrix_impl<T>&, future<int>&);

        inline void l(const matrix_impl<T>& a, const matrix_impl<T>& b, future<int>& ret){
            this->ctxt_select("1 from ambient as validation"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_pin(a); 
            this->block_2d_cycle_assign(b); 
        }
        
        inline void c(const matrix_impl<T>& a, const matrix_impl<T>& b, future<int>& ret){ // see paper for Reference Dongara 
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
                            ret.get_value() = 0; // test failed return 0 (bool false)
                        }
                    }
                }
            }
        }
    };

    // {{{ MKL LAPACK kernels

    template<typename T>
    struct svd : public kernel_unpinned< svd<T> > 
    {
        typedef void (svd::*F)(const matrix_impl<T>&, int&, int&, int&, matrix_impl<T>&, 
                          matrix_impl<T>&, matrix_impl<double>&);

        inline void l(const matrix_impl<T>& a, int& m, int& n, int& k, matrix_impl<T>& u, 
                      matrix_impl<T>& vt, matrix_impl<double>& s)
        {
            this->ctxt_select("1 from ambient as svd"); //if(!ctxt.involved()) return;
            this->block_outright_conditional_assign(s);
            this->block_2d_cycle_conditional_assign(a);
            this->block_2d_cycle_conditional_assign(u);
            this->block_2d_cycle_conditional_assign(vt);
        }

        inline void c(const matrix_impl<T>& a, int& m, int& n, int& k, matrix_impl<T>& u, 
                     matrix_impl<T>& vt, matrix_impl<double>& s)
        {
            // gs
            __A_TIME_C("ambient_svd_c_kernel"); 
        // Locals //
            int lda = ui_c_get_grid_dim(a).y*ui_c_get_mem_dim(a).y;
            int ldu = ui_c_get_grid_dim(u).y*ui_c_get_mem_dim(u).y;
            int ldvt = ui_c_get_grid_dim(vt).y*ui_c_get_mem_dim(vt).y;
            int info, lwork;
            T wkopt;
            T* work;
            double* rwork = new double[5*std::min(lda,ldu)]; // C - useless for double but need for complex 
            T* ad  = (T*)__a_solidify<T>(a);
            T* ud  = (T*)__a_solidify<T>(u);
            T* vtd = (T*)__a_solidify<T>(vt);
            double* sd = (double*)__a_solidify<T>(s);
            
        // Query and allocate the optimal workspace //
            lwork = -1; // C - Alex, netlib said -1 for the best workspace
            gesvd( "S", "S", &m, &n, ad, &lda, sd, ud, &ldu, vtd, &ldvt, &wkopt, &lwork, rwork, &info );
            lwork = OptimalSize(wkopt);
            work = (T*)malloc( lwork*sizeof(T) );
        // Compute SVD //
            gesvd( "S", "S", &m, &n, ad, &lda, sd, ud, &ldu, vtd, &ldvt, work, &lwork, rwork, &info );
        // Check for convergence //
            assert( info == 0 ); // otherwise the algorithm computing SVD failed to converge
            __a_disperse<T>(ud, u);
            __a_disperse<T>(vtd, vt);
            __a_disperse<T>(sd, s);
            free(work);
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct syev : public kernel_unpinned< syev<T> > 
    {
        typedef void (syev::*F)(matrix_impl<T>&, int&, matrix_impl<T>&);

        inline void l(matrix_impl<T>& a, int& m, matrix_impl<T>& w){
            this->ctxt_select("1 from ambient as syev"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_conditional_assign(a);
            this->block_2d_cycle_conditional_assign(w);
        }

        inline void c(matrix_impl<T>& a, int& m, matrix_impl<T>& w){
            __A_TIME_C("ambient_syev_c_kernel"); 
            int lda = ui_c_get_grid_dim(a).y*ui_c_get_mem_dim(a).y;
            int info, lwork = -1;
            double wkopt;
            double* work;
            double* ad = (double*)__a_solidify<T>(a);
            double* wd = (double*)__a_solidify<T>(w);
             
            dsyev_("V","U",&m,ad,&lda,wd,&wkopt,&lwork,&info);
            lwork = (int)wkopt;
            work = (double*)malloc( lwork*sizeof(double) );
            dsyev_("V","U",&m,ad,&lda,wd,work,&lwork,&info);
            assert( info == 0 ); // otherwise the algorithm computing SYEV failed to converge
            __a_disperse<T>(ad, a);
            __a_disperse<T>(wd, w);
            free(work); 
            __A_TIME_C_STOP
        }
    };
        
    template<typename T>
    struct heev : public kernel_unpinned< heev<T> > 
    {
        typedef void (heev::*F)(matrix_impl<T>&, const size_t&, matrix_impl<double>&);

        inline void l(matrix_impl<T>& a, const size_t& m, matrix_impl<double>& w){
            this->ctxt_select("1 from ambient as heev"); //if(!ctxt.involved()) return;
            this->block_2d_cycle_conditional_assign(a);
            this->block_2d_cycle_conditional_assign(w); // C - block_outright(w) is possible, if yes remove solidify and disperse for w 
        }

        inline void c(matrix_impl<T>& a, const size_t& m, matrix_impl<double>& w){
            // gs
            __A_TIME_C("ambient_heev_c_kernel"); 
            int lda = ui_c_get_grid_dim(a).y*ui_c_get_mem_dim(a).y;
            int info, lwork = -1;
        
            double wkopt;
            double* work;
            int am = (int)m; // for mkl (int*)
        
            double* ad = (double*)__a_solidify<double>(a);
            double* wd = (double*)__a_solidify<double>(w);
       
            dsyev_("V","U",&am,ad,&lda,wd,&wkopt,&lwork,&info);
            lwork = (int)wkopt;
            work = (double*)malloc( lwork*sizeof(double) );
            dsyev_("V","U",&am,ad,&lda,wd,work,&lwork,&info);
            assert(info == 0); // otherwise the algorithm computing HEEV failed to converge
            
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
         
            __a_disperse<double>(ad, a);
            __a_disperse<double>(wd, w);
            free(work);
            __A_TIME_C_STOP
        }
    };

    // }}}

} } }
#endif
