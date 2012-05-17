#ifndef __MAQUIS_TYPES_C_KERNELS_HPP__
#define __MAQUIS_TYPES_C_KERNELS_HPP__
// C - Note about the vectorization, the maximum of iteration is calculated outside the loop, to help at least the intel compiler to vectorize
#include <limits> // for the validation test
#include "utils/timings.h"
namespace ambient {

    #include "ambient/utils/numeric.h" // Blas/Lapack signature
    #include "ambient/utils/ceil.h"

    //#define AMBIENT_COMPUTATIONAL_TIMINGS
   
    #ifdef AMBIENT_COMPUTATIONAL_TIMINGS
        #define __A_TIME(name) static TimerPTH time(name); time.begin();
        #define __A_TIME_STOP time.end();
    #else
        #define __A_TIME(name) 
        #define __A_TIME_STOP 
    #endif

    template<typename T>
    inline size_t __a_get_limit_x(const T& a, size_t n = 0){
        if(n == 0) n = ui_c_get_dim(a).x;
        return std::min(ui_c_get_mem_dim(a).x, n-ctxt.get_block_id().x*ui_c_get_mem_dim(a).x);
    }

    template<ambient::models::v_model::revision&(*STATE)(const ambient::models::v_model::object&), typename T>
    inline size_t __a_get_limit_x(const T& a, size_t n = 0){
        if(n == 0) n = ui_c_get_dim(a).x;
        return std::min(ui_c_get_mem_dim<STATE>(a).x, n-ctxt.get_block_id().x*ui_c_get_mem_dim<STATE>(a).x);
    }

    template<typename T>
    inline size_t __a_get_limit_y(const T& a, size_t m = 0){
        if(m == 0) m = ui_c_get_dim(a).y;
        return std::min(ui_c_get_mem_dim(a).y, m-ctxt.get_block_id().y*ui_c_get_mem_dim(a).y);
    }

    template<ambient::models::v_model::revision&(*STATE)(const ambient::models::v_model::object&), typename T>
    inline size_t __a_get_limit_y(const T& a, size_t m = 0){
        if(m == 0) m = ui_c_get_dim<STATE>(a).y;
        return std::min(ui_c_get_mem_dim<STATE>(a).y, m-ctxt.get_block_id().y*ui_c_get_mem_dim<STATE>(a).y);
    }

    // C - w, alfa are not nested into dim2 because they are not 2D coordinates
    template <typename T>
    void __a_memcpy(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 const& dpos, maquis::types::p_dense_matrix_impl<T> const& src, T *sd, dim2 const& spos, size_t w, T alfa){
        memcpy(&dd[dpos.x*ui_c_get_mem_dim<ui_c_updated>(dest).y+dpos.y],
               &sd[spos.x*ui_c_get_mem_dim(src).y+spos.y],
               w*sizeof(T));
    }

    template <typename T>
    void __a_memscal(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 const& dpos, maquis::types::p_dense_matrix_impl<T> const& src, T *sd, dim2 const& spos, size_t w, T alfa){
        for(int z = 0; z < w; z++)
            dd[dpos.x*ui_c_get_mem_dim<ui_c_updated>(dest).y+dpos.y+z] += sd[spos.x*ui_c_get_mem_dim(src).y+spos.y + z]*alfa;
    }

    // V = double, S = size_t
    template<typename T>
    inline void __a_memptf(void (*ptf)(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 const& dpos, 
                                       maquis::types::p_dense_matrix_impl<T> const& src, T *sd, dim2 const& spos, 
                                       size_t w, T alfa),
                           maquis::types::p_dense_matrix_impl<T>& dest, dim2 dest_p, 
                           const maquis::types::p_dense_matrix_impl<T>& src, dim2 src_p, 
                           dim2 size, T alfa = 0.0)
    {
        __A_TIME("ambient_memptf_f_kernel");
        // the ouput (dest) must be a pinned p_dense_matrix

        T* dd = ui_c_updated(dest)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t di = ctxt.get_block_id().y * ui_c_get_mem_dim(dest).y;
        size_t dj = ctxt.get_block_id().x * ui_c_get_mem_dim(dest).x;
    
        if(ui_c_get_dim(dest).x - dest_p.x < size.x || ui_c_get_dim(dest).y - dest_p.y < size.y ||
           ui_c_get_dim(src).x - src_p.x   < size.x || ui_c_get_dim(src).y - src_p.y   < size.y) 
            maquis::cout << "Error: invalid memory movement" << std::endl;
    
        if( size.x == 0 || size.y == 0            ||
            di + ui_c_get_mem_dim(dest).y <= dest_p.y  || 
            dj + ui_c_get_mem_dim(dest).x <= dest_p.x  ||
            di >= dest_p.y + size.y               || 
            dj >= dest_p.x + size.x               ) return;

    // lets find dest-block starting point
        size_t offseti = (di < dest_p.y) ? dest_p.y % ui_c_get_mem_dim(dest).y : 0; di += offseti;
        size_t offsetj = (dj < dest_p.x) ? dest_p.x % ui_c_get_mem_dim(dest).x : 0; dj += offsetj;

        size_t si = di - dest_p.y + src_p.y;
        size_t sj = dj - dest_p.x + src_p.x;
    // lets find dest-block copy limits
        size_t sizey = __a_get_limit_y(dest, (dest_p.y + size.y)) - offseti;
        size_t sizex = __a_get_limit_x(dest, (dest_p.x + size.x)) - offsetj;

    // how many blocks do we need for this one
        int remaining = (int)sizey - (int)(ui_c_get_mem_dim(src).y - si % ui_c_get_mem_dim(src).y);
        size_t num_src_blocks = remaining > 0 ? __a_ceil( remaining / ui_c_get_mem_dim(src).y ) + 1 : 1;

        dim2 dpos,spos; size_t v, w;
        for(size_t j = 0; j < sizex; j++){
            w = sizey;
            dpos.y = offseti;
            dpos.x = offsetj + j;
            spos.y = si % ui_c_get_mem_dim(src).y;
            spos.x = (sj+j) % ui_c_get_mem_dim(src).x;
            for(int k = 0; k < num_src_blocks; k++){
                T* sd = ui_c_current(src)(si / ui_c_get_mem_dim(src).y + k, (sj+j) / ui_c_get_mem_dim(src).x);
                v = std::min(w, ui_c_get_mem_dim(src).y-spos.y);
                ptf(dest,dd,dpos,src,sd,spos,v,alfa);            
                w -= v;
                dpos.y += v;
                spos.y = 0;
            }
        }
        __A_TIME_STOP
    }

    template<typename T>
    inline void __a_memptf_reverse(void (*ptf)(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 const& dpos, 
                                               maquis::types::p_dense_matrix_impl<T> const& src, T *sd, dim2 const& spos, 
                                               size_t w, T alfa),
                                   maquis::types::p_dense_matrix_impl<T>& dest, dim2 dest_p, 
                                   const maquis::types::p_dense_matrix_impl<T>& src, dim2 src_p, 
                                   dim2 size, T alfa = 0.0)
    {
        __A_TIME("ambient_memptf_revers_f_kernel");
        // the input (src) must be a pinned p_dense_matrix

        T* sd = ui_c_current(src)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t si = ctxt.get_block_id().y * ui_c_get_mem_dim(src).y;
        size_t sj = ctxt.get_block_id().x * ui_c_get_mem_dim(src).x;
    
        if(ui_c_get_dim<ui_c_updated>(dest).x - dest_p.x < size.x || ui_c_get_dim<ui_c_updated>(dest).y - dest_p.y < size.y ||
           ui_c_get_dim(src).x - src_p.x   < size.x || ui_c_get_dim(src).y - src_p.y   < size.y) 
            maquis::cout << "Error: invalid memory movement" << std::endl;
    
        if( size.x == 0 || size.y == 0          ||
            si + ui_c_get_mem_dim(src).y <= src_p.y  || 
            sj + ui_c_get_mem_dim(src).x <= src_p.x  ||
            si >= src_p.y + size.y              || 
            sj >= src_p.x + size.x               ) return;

    // lets find src-block starting point
        size_t offseti = (si < src_p.y) ? src_p.y % ui_c_get_mem_dim(src).y : 0; si += offseti;
        size_t offsetj = (sj < src_p.x) ? src_p.x % ui_c_get_mem_dim(src).x : 0; sj += offsetj;

        size_t di = si - src_p.y + dest_p.y;
        size_t dj = sj - src_p.x + dest_p.x;
    // lets find src-block copy limits
        size_t sizey = __a_get_limit_y(src, (src_p.y + size.y)) - offseti;
        size_t sizex = __a_get_limit_x(src, (src_p.x + size.x)) - offsetj;

    // how many blocks do we need for this one
        int remaining = (int)sizey - (int)(ui_c_get_mem_dim<ui_c_updated>(dest).y - di % ui_c_get_mem_dim<ui_c_updated>(dest).y);
        size_t num_dest_blocks = remaining > 0 ? __a_ceil( remaining / ui_c_get_mem_dim<ui_c_updated>(dest).y ) + 1 : 1;

        dim2 dpos,spos; size_t v, w;
        for(size_t j = 0; j < sizex; j++){
            w = sizey;
            spos.y = offseti;
            spos.x = offsetj + j;
            dpos.y = di % ui_c_get_mem_dim<ui_c_updated>(dest).y;
            dpos.x = (dj+j) % ui_c_get_mem_dim<ui_c_updated>(dest).x;
            for(int k = 0; k < num_dest_blocks; k++){
                T* dd = ui_c_updated(dest)(di / ui_c_get_mem_dim<ui_c_updated>(dest).y + k, (dj+j) / ui_c_get_mem_dim<ui_c_updated>(dest).x);
                v = std::min(w, ui_c_get_mem_dim<ui_c_updated>(dest).y-dpos.y);
                ptf(dest,dd,dpos,src,sd,spos,v,alfa);            
                w -= v;
                spos.y += v;
                dpos.y = 0;
            }
        }
        __A_TIME_STOP
    }

    template<typename T>
    void print_pinned_block(T& a){
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
    
        double* ad = ui_c_current(a)(i, j);
        std::cout << " --(i,j)-- " << i << "," << j << std::endl;
        for(int ii = 0; ii < ui_c_get_mem_dim(a).y; ii++){
            for(int jj = 0; jj < ui_c_get_mem_dim(a).x; jj++){
                std::cout << ad[jj*ui_c_get_mem_dim(a).y + ii]<< " ";
            }     
            std::cout << " " << std::endl;
        }
        std::cout << " --------------------------- " << std::endl;
    }

    template<typename T>
    void gemm_inplace_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |  x  | 0 | 1 | 2 |  =  | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    //
    // partial reduce?..
    /////////////////////////////////////////////////////////////////////////
        int m   = ui_c_get_mem_dim(a).y;
        int n   = ui_c_get_mem_dim(b).x;
        int k   = ui_c_get_mem_dim(b).y;
        int lda = m;
        int ldb = k;
        int ldc = m;
        T alpha(1.0); 
        T beta(1.0);
    // a(i,j) => a(z,j) x b(j,i)  where z : [0,m)
    // current block of matrix a:
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
    // taking (j,i) of b:
        if(ui_c_get_grid_dim(b).y > j) while(i < ui_c_get_grid_dim(b).x){
            T* bd = ui_c_current(b)(j,i); // remote
    // multiplying with column of a:
            std::list<int> L;
            for(int z = 0; z < ui_c_get_grid_dim(a).y; z++) L.push_back(z);
            while(!L.empty()){
                std::list<int>::iterator zi = L.begin();
                while(zi != L.end()){
                    if(!ui_c_updated(a)(*zi,i).trylock()){ zi++; continue; }
                    T* ad = ui_c_current(a)(*zi,j);
                    T* cd = ui_c_updated(a)(*zi,i); // a(z,j) x b(j,i) => c(z,i)
                    gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
                    ui_c_updated(a)(*zi,i).unlock();
                    L.erase(zi++);
                }
            }
            i += ui_c_get_grid_dim(a).y;
        }
    }

    template<typename T>
    void gemm_one_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
        // gs // 0.2
        __A_TIME("ambient_gemm_one_c_kernel");
        (*(T*)ui_c_updated(c)(0,0)) = (*(T*)ui_c_current(a)(0,0))*(*(T*)ui_c_current(b)(0,0));
        __A_TIME_STOP
    }

    template<typename T>
    void gemm_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
        // gs // 0.2
        __A_TIME("ambient_gemm_c_kernel");
        if(ui_c_get_grid_dim(a) == 1 && ui_c_get_grid_dim(b) == 1){
            T* bd = ui_c_current(b)(0,0);
            T* ad = ui_c_current(a)(0,0);
            T* cd = ui_c_updated(c)(0,0);
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
            return;
        }else if(ui_c_get_mem_dim(a) != ui_c_get_mem_dim(b) || ui_c_get_mem_dim(a) != ui_c_get_mem_dim(c)){
            T* ad = (T*)models::solidify(a);
            T* bd = (T*)models::solidify(b);
            T* cd = (T*)models::solidify(c);
            int m   = ui_c_get_dim(a).y;
            int n   = ui_c_get_dim(b).x;
            int k   = ui_c_get_dim(b).y;
            int lda = ui_c_get_mem_dim(a).y*ui_c_get_grid_dim(a).y;
            int ldb = ui_c_get_mem_dim(b).y*ui_c_get_grid_dim(b).y;
            int ldc = ui_c_get_mem_dim(c).y*ui_c_get_grid_dim(c).y;
            T alpha(1.0); 
            T beta(1.0);
            gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
            models::disperse(cd, c);
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
    // a(i,j) => a(z,j) x b(j,i)  where z : [0,m)
    // current block of matrix a:
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;

        dim2 a_grid_dim = ui_c_get_grid_dim(a);
        dim2 b_grid_dim = ui_c_get_grid_dim(b);
    // taking (j,i) of b:
        if(b_grid_dim.y > j) while(i < b_grid_dim.x){
            T* bd = ui_c_current(b)(j,i); // remote
    // multiplying with column of a:
            std::list<int> L;
            for(int z = 0; z < a_grid_dim.y; z++) L.push_back(z);
            while(!L.empty()){
                std::list<int>::iterator zi = L.begin();
                while(zi != L.end()){
                    if(!ui_c_updated(c)(*zi,i).trylock()){ zi++; continue; }
                    T* ad = ui_c_current(a)(*zi,j);
                    T* cd = ui_c_updated(c)(*zi,i); // a(z,j) x b(j,i) => c(z,i)
                    gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
                    ui_c_updated(c)(*zi,i).unlock();
                    L.erase(zi++);
                }
            }
            i += a_grid_dim.y;
        }
        __A_TIME_STOP
    }

    template<typename T>
    void copy_c(maquis::types::p_dense_matrix_impl<T>& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a){
        // gs // 0.65
        __A_TIME("ambient_copy_c_kernel"); 
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        T* a_elements  = ui_c_current(a)(i,j);
        T* ac_elements = ui_c_updated(ac)(i,j);
        memcpy(ac_elements, a_elements, sizeof(T)*ui_c_get_mem_dim(a).y*ui_c_get_mem_dim(a).x);
        __A_TIME_STOP
    }

    template<typename T>
    void remove_rows_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& i_mark, const size_t& k){
        size_t numrows = ui_c_get_dim(a).y;
        size_t numcols = ui_c_get_dim(a).x;
        __a_memptf_reverse(&__a_memcpy<T>, a, dim2(0,0), a, dim2(0,0), dim2(numcols, i_mark));
        __a_memptf_reverse(&__a_memcpy<T>, a, dim2(0,i_mark), a, dim2(0,k+i_mark), dim2(numcols,numrows-k-i_mark));
    }

    template<typename T>
    void remove_cols_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& j_mark, const size_t& k){
        size_t numrows = ui_c_get_dim(a).y;
        size_t numcols = ui_c_get_dim(a).x;
        __a_memptf_reverse(&__a_memcpy<T>, a, dim2(0,0), a, dim2(0,0), dim2(j_mark, numrows));
        __a_memptf_reverse(&__a_memcpy<T>, a, dim2(j_mark,0), a, dim2(k+j_mark,0), dim2(numcols-k-j_mark,numrows));
    }

    template<typename T>
    void resize_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& om, const size_t& on){
        __a_memptf_reverse(&__a_memcpy<T>, a, dim2(0,0), a, dim2(0,0), dim2(std::min(n,on), std::min(m,om)));
    }

    template<typename T>
    void sqrt_diagonal_c(pinned maquis::types::p_dense_matrix_impl<T>& a){
        T* ad = ui_c_current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        T* sd = ui_c_updated(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t size = ui_c_get_mem_dim(a).y;
        for(int i=0; i < size; i++)
            sd[i] = sqrt(ad[i]);
    }

    template<typename T>
    void exp_diagonal_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const T& alfa){
        T* ad = ui_c_current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        T* sd = ui_c_updated(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t size = ui_c_get_mem_dim(a).y;
        for(int i=0; i < size; i++)
            sd[i] = exp(ad[i]*alfa);
    }

    template<typename T>
    void exp_diagonal_rc_c(maquis::types::p_dense_matrix_impl< std::complex<T> >& e, pinned const maquis::types::p_dense_matrix_impl<T>& a, const std::complex<T>& alfa){
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        std::complex<T>* ed = ui_c_updated(e)(i, j);
        T* ad = ui_c_current(a)(i, j);
        size_t size = ui_c_get_mem_dim(e).y;
        for(int i=0; i < size; i++)
            ed[i] = exp(ad[i]*alfa);
    }

    template<typename T>
    void push_back_sqr_gt_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, std::vector<T>*& ac){
        // gs
        __A_TIME("ambient_push_back_sqr_gt_c_kernel"); 
        double* ad = ui_c_current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t sizey = __a_get_limit_y(a);
        for(int i=0; i < sizey; i++){
            double v = std::abs(ad[i]);
            if(v > 1e-10) ac->push_back(v*v);
        }
        __A_TIME_STOP
    }

    template<typename T>
    void cast_to_dense_c(std::vector<T>*& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
        // gs
        __A_TIME("ambient_cast_to_dense_c_kernel"); 
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        int yi = ui_c_get_mem_dim(a).y*i; // conversion cartersia coordinates dense / p_dense
        int xj = ui_c_get_mem_dim(a).x*j; 
        size_t offset;
        T* ad = ui_c_current(a)(i,j);
       
        size_t sizex = __a_get_limit_x(a, n);
        size_t sizey = __a_get_limit_y(a, m);
        size_t lda = ui_c_get_mem_dim(a).y;
    
        for(int jj=0; jj < sizex; ++jj){
            offset = yi + (xj+jj)*m;
            memcpy((void*)&(*ac)[offset],(void*)&ad[jj*lda], sizey*sizeof(T));  
        }
        __A_TIME_STOP
    }

    template <typename T>
    void cast_to_p_dense_c(const std::vector<T>*& ac, pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& lda){
        int offset,size_y(ui_c_get_mem_dim(a).y),size_x(ui_c_get_mem_dim(a).x);
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        int yi = ui_c_get_mem_dim(a).y*i; // conversion cartersia coordinates dense / p_dense
        int xj = ui_c_get_mem_dim(a).x*j; 
        T* ad = ui_c_updated(a)(i,j);
    
        //operator ?, case 1 matrix is lower than one work group, case 2 several work groups or fit in x direction
        if(j+1 == ui_c_get_grid_dim(a).x)
            size_x = (ui_c_get_mem_dim(a).x > n) ? n : (n - (ui_c_get_grid_dim(a).x-1)*ui_c_get_mem_dim(a).x);
    
        for(int ii=0; ii < size_x; ++ii){
            offset = yi + (xj+ii)*lda; //lda because possible resize;
            if(i+1 == ui_c_get_grid_dim(a).y)
               size_y = (ui_c_get_mem_dim(a).y > m) ? m : (m - (ui_c_get_grid_dim(a).y-1)*ui_c_get_mem_dim(a).y);
            memcpy((void*)&ad[ii*ui_c_get_mem_dim(a).x],(void*)&(*ac)[offset], size_y*sizeof(T)); // y direction 
        }
    }

    template <typename T>
    void reshape_l2r_c(const maquis::types::p_dense_matrix_impl<T>& left, pinned maquis::types::p_dense_matrix_impl<T>& right,
                       const size_t& left_offset, const size_t& right_offset, 
                       const size_t& sdim, const size_t& ldim, const size_t& rdim)
    { // gs
        __A_TIME("ambient_reshape_l2r_c_kernel"); 
        __a_memptf(&__a_memcpy<T>, right, dim2(0,0), right, dim2(0,0), dim2(ui_c_get_dim(right).x,ui_c_get_dim(right).y)); // refreshing updated memory
        for(size_t ss = 0; ss < sdim; ++ss){
            __a_memptf(&__a_memcpy<T>, right, dim2(ss*rdim + right_offset, 0), 
                       left,  dim2(0, ss*ldim + left_offset), 
                       dim2( rdim, ldim ));
        }
        __A_TIME_STOP
    }

    template <typename T>
    void reshape_r2l_c(pinned maquis::types::p_dense_matrix_impl<T>& left, const maquis::types::p_dense_matrix_impl<T>& right,
                       const size_t& left_offset, const size_t& right_offset, 
                       const size_t& sdim, const size_t& ldim, const size_t& rdim)
    { // gs
        __A_TIME("ambient_reshape_r2l_c_kernel"); 
        __a_memptf(&__a_memcpy<T>, left, dim2(0,0), left, dim2(0,0), dim2(ui_c_get_dim(left).x,ui_c_get_dim(left).y)); // refreshing updated memory
        for(size_t ss = 0; ss < sdim; ++ss)
            __a_memptf(&__a_memcpy<T>, left,  dim2(0, ss*ldim + left_offset), 
                       right, dim2(ss*rdim + right_offset,0), 
                       dim2( rdim, ldim ));
        __A_TIME_STOP
    }

    template <typename T>
    void rb_tensor_mpo_c(pinned maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                         const size_t& out_offset, const size_t& in_offset, 
                         const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
    { // gs
        __A_TIME("ambient_rb_tensor_mpo_c_kernel"); 
        __a_memptf(&__a_memcpy<T>, out, dim2(0,0), out, dim2(0,0), dim2(ui_c_get_dim(out).x,ui_c_get_dim(out).y)); // refreshing updated memory
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
            for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                T* alfad = ui_c_current(alfa)(ss1/ui_c_get_mem_dim(alfa).y, ss2/ui_c_get_mem_dim(alfa).x);
                T  alfa_t = alfad[ss1%ui_c_get_mem_dim(alfa).y + ui_c_get_mem_dim(alfa).y*(ss2%ui_c_get_mem_dim(alfa).x)];
                __a_memptf(&__a_memscal<T>, out, dim2(out_offset + ss2*rdim, 0),
                                            in,  dim2(in_offset + ss1*rdim, 0),
                                            dim2(rdim, ldim), alfa_t);
            }
        __A_TIME_STOP
    }

    template <typename T>
    void lb_tensor_mpo_c(pinned maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                         const size_t& out_offset, const size_t& in_offset, 
                         const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
    { // gs
        __A_TIME("ambient_lb_tensor_mpo_c_kernel"); 
        __a_memptf(&__a_memcpy<T>, out, dim2(0,0), out, dim2(0,0), dim2(ui_c_get_dim(out).x,ui_c_get_dim(out).y)); // refreshing updated memory
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
            for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                T* alfad = ui_c_current(alfa)(ss1/ui_c_get_mem_dim(alfa).y, ss2/ui_c_get_mem_dim(alfa).x);
                T  alfa_t = alfad[ss1%ui_c_get_mem_dim(alfa).y + ui_c_get_mem_dim(alfa).y*(ss2%ui_c_get_mem_dim(alfa).x)];
                __a_memptf(&__a_memscal<T>, out, dim2(0, out_offset + ss2*ldim),
                                            in,  dim2(0, in_offset + ss1*ldim),
                                            dim2(rdim, ldim), alfa_t);
            }
        __A_TIME_STOP
    }

    template<typename T>
    void scalar_norm_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, T*& norm){
        // gs
        __A_TIME("ambient_scalar_norm_c_kernel"); 
        T summ = 0;
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        T* ad = ui_c_current(a)(i,j);
        size_t lda = ui_c_get_mem_dim(a).y;
        size_t sizey = __a_get_limit_y(a, m);
        size_t sizex = __a_get_limit_x(a, n);
        for(size_t ii=0; ii < sizey; ii++)
            for(size_t jj=0; jj < sizex; jj++)
                summ += ad[ii+jj*lda]*ad[ii+jj*lda];
    
        *norm += summ;
        __A_TIME_STOP
    }

    template<typename T>
    void scalar_overlap_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, const size_t& m, const size_t& n, T*& overlap){
        // gs
        __A_TIME("ambient_scalar_overlap_c_kernel"); 
        T summ = 0;
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        T* ad = ui_c_current(a)(i,j);
        T* bd = ui_c_current(b)(i,j);
        size_t lda = ui_c_get_mem_dim(a).y;
        size_t sizey = __a_get_limit_y(a);
        size_t sizex = __a_get_limit_x(a);
        for(size_t ii=0; ii < sizey; ii++)
            for(size_t jj=0; jj < sizex; jj++)
                summ += ad[ii+jj*lda]*bd[ii+jj*lda];
        *overlap += summ;
        __A_TIME_STOP
    }

    template<typename T>
    void add_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
        // gs
        __A_TIME("ambient_add_c_kernel"); 
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        T* ad = ui_c_current(a)(i, j);
        T* bd = ui_c_current(b)(i, j);
        T* ar = ui_c_updated(a)(i, j);
        size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
        for(size_t k = 0; k < size; k++)
            ar[k] = ad[k] + bd[k];
        __A_TIME_STOP
    }

    template<typename T>
    void sub_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
        // gs
        __A_TIME("ambient_sub_c_kernel"); 
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        double* ad = ui_c_current(a)(i, j);
        double* bd = ui_c_current(b)(i, j);
        double* ar = ui_c_updated(a)(i, j);
        size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
        for(size_t k = 0; k < size; k++)
            ar[k] = ad[k] + (-1)*bd[k];
        __A_TIME_STOP
    }

    template<typename T>
    void scale_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const double*& t){
        // gs
        __A_TIME("ambient_scale_c_kernel"); 
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        T* ad = ui_c_current(a)(i, j);
        T* ar = ui_c_updated(a)(i, j);

        size_t sizey = __a_get_limit_y(a, m);
        size_t sizex = __a_get_limit_x(a, n);
        size_t lda = ui_c_get_mem_dim(a).y;
        for(size_t jj=0; jj < sizex; jj++)
        for(size_t ii=0; ii < sizey; ii++)
        ar[jj*lda+ii] = ad[jj*lda+ii] * (*t);
        __A_TIME_STOP
    }

    template<typename T, typename D>
    void gemm_diagonal_lhs_c(const maquis::types::p_dense_matrix_impl<D>& a_diag, pinned const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c,
            const size_t& m, const size_t& n, const size_t& k){
        // gs
        __A_TIME("ambient_gemm_diagonal_lhs_c_kernel"); 
        size_t sizey = __a_get_limit_y(a_diag, m);
        int j = ctxt.get_block_id().y*ui_c_get_mem_dim(b).y;
        int size = ui_c_get_mem_dim(b).x;
        int lda  = sizeof(T)/sizeof(D)*ui_c_get_mem_dim(b).y;
        int ONE  = 1;
        D* bd = ui_c_current(b)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        D* cd = ui_c_updated(c)(ctxt.get_block_id().y, ctxt.get_block_id().x);
    
        for(int jj = 0 ; jj < sizey; jj++){
             D* alpha = ui_c_current(a_diag)((j+jj)/ui_c_get_mem_dim(a_diag).y,0);
    	     axpy(&size, &alpha[(j+jj)%ui_c_get_mem_dim(a_diag).y], &bd[jj], &lda, &cd[jj], &lda);
    	     if(sizeof(T) != sizeof(D)) axpy(&size, &alpha[(j+jj)%ui_c_get_mem_dim(a_diag).y], &bd[jj+1], &lda, &cd[jj], &lda); // for complex
        }
        __A_TIME_STOP
    }

    template<typename T, typename D>
    void gemm_diagonal_rhs_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<D>& b_diag, maquis::types::p_dense_matrix_impl<T>& c,
            const size_t& m, const size_t& n, const size_t& k){
        // gs
        __A_TIME("ambient_gemm_diagonal_rhs_c_kernel"); 
        size_t sizex = std::min((ctxt.get_block_id().x+1)*ui_c_get_mem_dim(b_diag).y, n)-ctxt.get_block_id().x*ui_c_get_mem_dim(b_diag).y;

        int j = ctxt.get_block_id().x*ui_c_get_mem_dim(a).x;
        int size = sizeof(T)/sizeof(D)*ui_c_get_mem_dim(a).y; // for the case of complex
        int ONE = 1;
        D* ad = ui_c_current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        D* cd = ui_c_updated(c)(ctxt.get_block_id().y, ctxt.get_block_id().x);
    
        for(int jj = 0 ; jj < sizex; jj++){
    	    D* alpha = ui_c_current(b_diag)((j+jj)/ui_c_get_mem_dim(b_diag).y,0);
    	    axpy(&size, &alpha[(j+jj)%ui_c_get_mem_dim(b_diag).y], &ad[jj*ui_c_get_mem_dim(a).y], &ONE, &cd[jj*ui_c_get_mem_dim(c).y], &ONE);
        }
        __A_TIME_STOP
    }

    template<typename T>
    void trace_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& n, T*& trace){
        // gs
        __A_TIME("ambient_trace_c_kernel"); 
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        size_t ld = ui_c_get_mem_dim(a).y;
        size_t sd = ui_c_get_mem_dim(a).x;
        T* ad = ui_c_current(a)(i,j);
    
        if((i+1)*ld <= j*sd) return;
        if(i*ld >= (j+1)*sd) return;
        size_t sizex = std::min(n,(j+1)*sd);
        for(size_t jj = j*sd; jj < sizex; jj++){
            if(i*ld > jj) continue;
            if((i+1)*ld <= jj) continue;
           *trace += ad[jj % ld + (jj%sd)*ld];
        }
        __A_TIME_STOP
    }

    template<typename T>
    void transpose_c(pinned maquis::types::p_dense_matrix_impl<T>& m){ // we need to reset dims also for non-square matrices
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        T* td = ui_c_updated(m)(i,j);
        T* od = ui_c_current(m)(j,i);
    
        for(size_t i = 0; i < ui_c_get_mem_dim(m).y; ++i){
            for(size_t j=0; j < ui_c_get_mem_dim(m).x; ++j){
                td[j+i*ui_c_get_mem_dim(m).y] = od[i+j*ui_c_get_mem_dim(m).y];
            }
        }
    }

    template<typename T>
    void transpose_out_c(pinned const maquis::types::p_dense_matrix_impl<T>& m, maquis::types::p_dense_matrix_impl<T>& t){
        // gs
        __A_TIME("ambient_transpose_out_c_kernel"); 
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        T* od = ui_c_current(m)(i,j);
        T* td = ui_c_updated(t)(j,i);
   
        size_t sizey = ui_c_get_mem_dim(m).y;
        size_t sizex = ui_c_get_mem_dim(m).x;
        for(size_t i = 0; i < sizey; ++i){
            for(size_t j=0; j < sizex; ++j){
                td[j+i*sizey] = od[i+j*sizey];
            }
        }
        __A_TIME_STOP
    }

    template<typename T>
    void init_value_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const T& value){
        __A_TIME("ambient_init_value_c_kernel"); 
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;

        T* ad = ui_c_updated(a)(i,j);
        size_t sizey = __a_get_limit_y(a, m);
        size_t sizex = __a_get_limit_x(a, n);
    
        for(size_t j=0; j < sizex; ++j){
            for(size_t i = 0; i < sizey; ++i){
                ad[i+j*ui_c_get_mem_dim(a).y] = value; // not a memset due to complex
            }
        }
        __A_TIME_STOP
    }

    template<typename T> void randomize(T* ad){ *ad = drand48(); }
    template<typename T> void randomize(std::complex<T>* ad){
        (*ad).real() = drand48();
        (*ad).imag() = drand48();
    }

    template<typename T>
    void init_random_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
        __A_TIME("ambient_init_random_c_kernel"); 
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        size_t ld = ui_c_get_mem_dim(a).y;
        size_t sd = ui_c_get_mem_dim(a).x;
        T* ad = ui_c_updated(a)(i,j);
        size_t sizey = __a_get_limit_y(a, m);
        size_t sizex = __a_get_limit_x(a, n);
       
        for(size_t jj = 0; jj < sizex; jj++){
            for(size_t ii = 0; ii < sizey; ii++)
                randomize((ad+(jj*ld+ii)));
        }
        __A_TIME_STOP
    }

    template<typename T>
    void init_identity_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
        __A_TIME("ambient_init_identity_c_kernel"); 
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        size_t ld = ui_c_get_mem_dim(a).y;
        size_t sd = ui_c_get_mem_dim(a).x;
        T* ad = ui_c_updated(a)(i,j);
        if((i+1)*ld <= j*sd) return;
        if(i*ld >= (j+1)*sd) return;
        size_t sizex = std::min(m,n); // respecting borders
        sizex = std::min(sizex,(j+1)*sd);
        for(size_t jj = j*sd; jj < sizex; jj++){
            if(i*ld > jj) continue;
            if((i+1)*ld <= jj) continue;
            ad[jj % ld + (jj%sd)*ld] = 1.;
        }
        __A_TIME_STOP
    }

    template <typename T> 
    void validation_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, int*& ret){ // see paper for Reference Dongara 
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        T* ad = ui_c_current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x); 
        T* bd = ui_c_current(b)(ctxt.get_block_id().y, ctxt.get_block_id().x); 
        double res(0.0); 
        double epsilon = std::numeric_limits<double>::epsilon();
        size_t position_x(0),position_y(0),position_xy(0); 
   
        for(size_t ii=0; ii < ui_c_get_mem_dim(a).y; ++ii){
            for(size_t jj=0; jj < ui_c_get_mem_dim(a).x; ++jj){
                position_x = j*ui_c_get_mem_dim(a).x+jj;
                position_y = i*ui_c_get_mem_dim(a).y+ii;
                if(position_x < ui_c_get_dim(a).x && position_y < ui_c_get_dim(a).y){
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

    // {{{ MKL LAPACK kernels

    template<typename T>
    void svd_c(const maquis::types::p_dense_matrix_impl<T>& a, int& m, int& n, maquis::types::p_dense_matrix_impl<T>& u, maquis::types::p_dense_matrix_impl<T>& vt, maquis::types::p_dense_matrix_impl<double>& s){
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
            ad  = (T*)models::solidify(a);
            ud  = (T*)models::solidify(u);
            vtd = (T*)models::solidify(vt);
            sd  = (double*)models::solidify(s);
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
            models::disperse(ud, u);
            models::disperse(vtd, vt);
            models::disperse(sd, s);
        }
        free(work);
        __A_TIME_STOP
    }

    template<typename T>
    void syev_c(maquis::types::p_dense_matrix_impl<T>& a, int& m, maquis::types::p_dense_matrix_impl<T>& w){
        int lda = ui_c_get_grid_dim(a).y*ui_c_get_mem_dim(a).y;
        int info, lwork = -1;
        double wkopt;
        double* work;
        double* ad = (double*)models::solidify(a);
        double* wd = ui_c_get_grid_dim(a) == 1 ? ui_c_updated(w)(0,0) : (double*)models::solidify(w);
         
        dsyev_("V","U",&m,ad,&lda,wd,&wkopt,&lwork,&info);
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        dsyev_("V","U",&m,ad,&lda,wd,work,&lwork,&info);
    
        if( info > 0 ) {
            printf( "The algorithm computing SYEV failed to converge.\n" );
            exit( 1 );
        }
    
        models::disperse(ad, a);
        if(ui_c_get_grid_dim(a) != 1) 
            models::disperse(wd, w);
        free(work); 
    }

    template<typename T>
    void heev_c(maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, maquis::types::p_dense_matrix_impl<T>& w){
        // gs
        __A_TIME("ambient_heev_c_kernel"); 
        int lda = ui_c_get_grid_dim(a).y*ui_c_get_mem_dim(a).y;
        int info, lwork = -1;
    
        double wkopt;
        double* work;
        int am = (int)m; // for mkl (int*)

        double* ad = (double*)models::solidify(a);
        double* wd = ui_c_get_grid_dim(a) == 1 ? ui_c_updated(w)(0,0) : (double*)models::solidify(w);
    
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
     
        models::disperse(ad, a);
        if(ui_c_get_grid_dim(a) != 1) 
            models::disperse(wd, w);
        free(work);
        __A_TIME_STOP
    }

    // }}}

    // {{{ strassen multiplication supplementary kernels

    template<typename T>
    void gemm_strassen_gad_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& ai, const size_t& aj, 
                             maquis::types::p_dense_matrix_impl<T>& r, const size_t& n)
    {
        size_t bi, bj;
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;

        size_t xi = i - (int)(ai / ui_c_get_mem_dim(a).y);
        size_t xj = j - (int)(aj / ui_c_get_mem_dim(a).x);

        int qr_grid_dim = (int)(n / ui_c_get_mem_dim(a).x)/2;
        dim2 qr = dim2((int)(xj / qr_grid_dim),(int)(xi / qr_grid_dim));

        // a11 + a12,  a12 - a22
        // a21 - a11,  a22 + a21
        if(qr.y == 0 && qr.x == 0){
            bi = 0 * qr_grid_dim + xi % qr_grid_dim;
            bj = 1 * qr_grid_dim + xj % qr_grid_dim;
        }else if(qr.y == 0 && qr.x == 1){
            bi = 1 * qr_grid_dim + xi % qr_grid_dim;
            bj = 1 * qr_grid_dim + xj % qr_grid_dim;
        }else if(qr.y == 1 && qr.x == 0){
            bi = 0 * qr_grid_dim + xi % qr_grid_dim;
            bj = 0 * qr_grid_dim + xj % qr_grid_dim;
        }else if(qr.y == 1 && qr.x == 1){
            bi = 1 * qr_grid_dim + xi % qr_grid_dim;
            bj = 0 * qr_grid_dim + xj % qr_grid_dim;
        }

        T* ad = ui_c_current(a)(i , j );
        T* rr = ui_c_updated(r)(xi, xj);
        T* bd = ui_c_current(a)(bi, bj);

        size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
        if(qr.x == qr.y) for(size_t k = 0; k < size; k++)
            rr[k] = ad[k] + bd[k];
        else for(size_t k = 0; k < size; k++)
            rr[k] = ad[k] - bd[k];
    }

    template<typename T>
    void gemm_strassen_dad_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& ai, const size_t& aj, 
                                    const maquis::types::p_dense_matrix_impl<T>& b, const size_t& bi, const size_t& bj, 
                                          maquis::types::p_dense_matrix_impl<T>& r, const size_t& n)
    { // a11 + a22,  b11 + b22
        int qr_grid_dim = (int)(n / ui_c_get_mem_dim(a).x)/2;
        int ci = ctxt.get_block_id().y;
        int cj = ctxt.get_block_id().x;
        int i = ci - (int)(ai / ui_c_get_mem_dim(a).y);
        int j = cj - (int)(aj / ui_c_get_mem_dim(a).x);

        if((int)(i / qr_grid_dim) != 0 || 
           (int)(j / qr_grid_dim) != 0) 
            return; // quick exit

        T* ad  = ui_c_current(a)(ci, cj);
        T* add = ui_c_current(a)(ci + qr_grid_dim, 
                            cj + qr_grid_dim);
        T* bd  = ui_c_current(b)( i + (int)(bi / ui_c_get_mem_dim(b).y), 
                             j + (int)(bj / ui_c_get_mem_dim(b).x));
        T* bdd = ui_c_current(b)( i + (int)(bi / ui_c_get_mem_dim(b).y) + qr_grid_dim,
                             j + (int)(bj / ui_c_get_mem_dim(b).x) + qr_grid_dim);
        T* arr = ui_c_updated(r)(i, j);
        T* brr = ui_c_updated(r)(i, j + qr_grid_dim);

        size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
        for(size_t k = 0; k < size; k++) arr[k] = ad[k] + add[k];
        for(size_t k = 0; k < size; k++) brr[k] = bd[k] + bdd[k];
    }

    template<typename T>
    void add_sum_submx_c(const  maquis::types::p_dense_matrix_impl<T>& a, const size_t& ai, const size_t& aj, 
                         const  maquis::types::p_dense_matrix_impl<T>& b, const size_t& bi, const size_t& bj, 
                         pinned maquis::types::p_dense_matrix_impl<T>& c, const size_t& ci, const size_t& cj, 
                         const size_t& n)
    { // c +=  a + b
        int i = ctxt.get_block_id().y - (int)(ci / ui_c_get_mem_dim(c).y);
        int j = ctxt.get_block_id().x - (int)(cj / ui_c_get_mem_dim(c).x);

        T* ad  = ui_c_current(a)( i + (int)(ai / ui_c_get_mem_dim(a).y), 
                             j + (int)(aj / ui_c_get_mem_dim(a).x));
        T* bd  = ui_c_current(b)( i + (int)(bi / ui_c_get_mem_dim(b).y), 
                             j + (int)(bj / ui_c_get_mem_dim(b).x));
        T* cd  = ui_c_current(c)( i + (int)(ci / ui_c_get_mem_dim(c).y), 
                             j + (int)(cj / ui_c_get_mem_dim(c).x));

        size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
        for(size_t k = 0; k < size; k++) cd[k] += ad[k] + bd[k];
    }

    template<typename T>
    void add_dif_submx_c(const  maquis::types::p_dense_matrix_impl<T>& a, const size_t& ai, const size_t& aj, 
                         const  maquis::types::p_dense_matrix_impl<T>& b, const size_t& bi, const size_t& bj, 
                         pinned maquis::types::p_dense_matrix_impl<T>& c, const size_t& ci, const size_t& cj, 
                         const size_t& n)
    { // c +=  a - b
        int i = ctxt.get_block_id().y - (int)(ci / ui_c_get_mem_dim(c).y);
        int j = ctxt.get_block_id().x - (int)(cj / ui_c_get_mem_dim(c).x);

        T* ad  = ui_c_current(a)( i + (int)(ai / ui_c_get_mem_dim(a).y), 
                             j + (int)(aj / ui_c_get_mem_dim(a).x));
        T* bd  = ui_c_current(b)( i + (int)(bi / ui_c_get_mem_dim(b).y), 
                             j + (int)(bj / ui_c_get_mem_dim(b).x));
        T* cd  = ui_c_current(c)( i + (int)(ci / ui_c_get_mem_dim(c).y), 
                             j + (int)(cj / ui_c_get_mem_dim(c).x));

        size_t size = ui_c_get_mem_dim(a).x*ui_c_get_mem_dim(a).y;
        for(size_t k = 0; k < size; k++) cd[k] += ad[k] - bd[k];
    }

    template<typename T>
    void gemm_submx_c(pinned const  maquis::types::p_dense_matrix_impl<T>& a, const size_t& ai, const size_t& aj, 
                             const  maquis::types::p_dense_matrix_impl<T>& b, const size_t& bi, const size_t& bj, 
                                    maquis::types::p_dense_matrix_impl<T>& c, const size_t& ci, const size_t& cj, 
                                    const size_t& size)
    { // atomic 1 block gemm
        T alpha, beta;
        int m, n, k, lda, ldb, ldc;
        m = n = k = lda = ldb = ldc = size;
        alpha = beta = 1.0; 

        size_t i = ctxt.get_block_id().y - (int)(ai / ui_c_get_mem_dim(a).y);
        size_t j = ctxt.get_block_id().x - (int)(aj / ui_c_get_mem_dim(a).x);

        T* ad  = ui_c_current(a)( i + (int)(ai / ui_c_get_mem_dim(a).y), 
                             j + (int)(aj / ui_c_get_mem_dim(a).x));
        T* bd  = ui_c_current(b)( i + (int)(bi / ui_c_get_mem_dim(b).y), 
                             j + (int)(bj / ui_c_get_mem_dim(b).x));
        T* cd  = ui_c_updated(c)( i + (int)(ci / ui_c_get_mem_dim(c).y), 
                             j + (int)(cj / ui_c_get_mem_dim(c).x));
        gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
    }

    // }}}

}
#endif
