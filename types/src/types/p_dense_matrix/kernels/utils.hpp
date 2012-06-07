#ifndef __MAQUIS_TYPES_KERNELS_UTILS_HPP__
#define __MAQUIS_TYPES_KERNELS_UTILS_HPP__

#include <limits>
#include "utils/timings.h"

extern "C" {
    double ddot_(const int*, const double*, const int*, const double*, const int*);
}

namespace ambient {

    #include "ambient/utils/numeric.h" // BLAS/LAPACK prototypes
    #include "ambient/utils/ceil.h"

    //#define AMBIENT_COMPUTATIONAL_TIMINGS
    //#define AMBIENT_CHECK_BOUNDARIES
   
    #ifdef AMBIENT_COMPUTATIONAL_TIMINGS
        #define __A_TIME(name) static TimerPTH time(name); time.begin();
        #define __A_TIME_STOP time.end();
    #else
        #define __A_TIME(name) 
        #define __A_TIME_STOP 
    #endif

    // {{{ continuous memory mangling
    template<typename V, typename T>
    inline void* __a_solidify(const iteratable<T>& o){
        using ambient::models::velvet::layout;

        fast_revision& r = o.ui_c_revision_0();
        layout& l = r.get_layout();
        size_t iterator = 0;
        char* memory = NULL;
        size_t stride = l.get_mem_dim().y*sizeof(V);
        size_t block = l.get_mem_dim().x;
        dim2 grid = l.get_grid_dim();

        for(size_t x=0; x < grid.x; x++)
            for(size_t xx=0; xx < block; xx++)
                for(size_t y=0; y < grid.y; y++){
                    memory = (char*)realloc(memory, (iterator+1)*stride);
                    memcpy(memory+iterator*stride,                         // copy to
                           &((char*)r(x,y))[xx*stride],                    // copy from
                           stride);                                        // of size
                    iterator++;
                }
        return memory;
    }

    template<typename V, typename T>
    inline void __a_disperse(void* data, iteratable<T>& o){
        using ambient::models::velvet::layout;
        using ambient::models::velvet::history;

        slow_revision& r = o.ui_c_revision_1();
        layout& l = r.get_layout();
        char* memory = (char*)data;
        size_t iterator = 0;
        size_t stride = l.get_mem_dim().y*sizeof(V);
        size_t block = l.get_mem_dim().x;
        dim2 grid = l.get_grid_dim();

        for(size_t x=0; x < grid.x; x++)
            for(size_t xx=0; xx < block; xx++)
                for(size_t y=0; y < grid.y; y++){
                    memcpy(&((char*)r(x,y))[xx*stride],                    // copy to
                           memory,                                         // copy from
                           stride);                                        // of size
                    memory += stride;
                }
        free(data);
    }
    // }}}

    template<typename T>
    inline size_t __a_get_limit_x(const T& a, size_t n = 0){
        if(n == 0) n = ui_c_get_dim(a).x;
        return std::min(ui_c_get_mem_dim(a).x, n-ctxt.get_block_id().x*ui_c_get_mem_dim(a).x);
    }

    template<typename T>
    inline size_t __a_get_limit_y(const T& a, size_t m = 0){
        if(m == 0) m = ui_c_get_dim(a).y;
        return std::min(ui_c_get_mem_dim(a).y, m-ctxt.get_block_id().y*ui_c_get_mem_dim(a).y);
    }

    template <typename T> inline T __a_dot(T* a, T* b, int size){
        T summ(0);
        for(size_t k=0; k < size; k++)
           summ += a[k]*b[k];
        return summ;
    }

    inline double __a_dot(double* a, double* b, int size){
        int ONE = 1;
        return ddot_(&size, a, &ONE, b, &ONE);
    }

    template<typename T>
    inline void __a_copy(T* dst, T* src, int size){
        memcpy(dst, src, size*sizeof(T));
    }

    template <typename T>
    inline void __a_memcpy(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 dpos, const maquis::types::p_dense_matrix_impl<T>& src, T *sd, dim2 spos, size_t w, T alfa){
        __a_copy(&dd[dpos.x*ui_c_get_mem_dim_u(dest).y+dpos.y],
                 &sd[spos.x*ui_c_get_mem_dim(src).y+spos.y],
                 w);
    }

    template <typename T>
    inline void __a_memscal(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 dpos, const maquis::types::p_dense_matrix_impl<T>& src, T *sd, dim2 spos, size_t w, T alfa){
        for(int z = 0; z < w; z++)
            dd[dpos.x*ui_c_get_mem_dim_u(dest).y+dpos.y + z] += sd[spos.x*ui_c_get_mem_dim(src).y+spos.y + z]*alfa;
    }

    template<typename T, void(*PTF)(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 dpos,
                                    const maquis::types::p_dense_matrix_impl<T>& src, T *sd, dim2 spos, 
                                    size_t w, T alfa)>
    inline void __a_memptf_atomic(maquis::types::p_dense_matrix_impl<T>& dest, dim2 dest_p, 
                                  const maquis::types::p_dense_matrix_impl<T>& src, dim2 src_p, 
                                  dim2 size, T alfa = 0.0)
    {
        __A_TIME("ambient_memptf_f_atomic_kernel");
        // the ouput (dest) must be a pinned p_dense_matrix

#ifdef AMBIENT_CHECK_BOUNDARIES
        if(ui_c_get_dim(dest).x - dest_p.x < size.x || ui_c_get_dim(dest).y - dest_p.y < size.y ||
           ui_c_get_dim(src).x - src_p.x   < size.x || ui_c_get_dim(src).y - src_p.y   < size.y) 
            maquis::cout << "Error: invalid memory movement" << std::endl;
#endif
    
        if( size.x == 0 || size.y == 0            ||
            ui_c_get_mem_dim(dest).y <= dest_p.y  || 
            ui_c_get_mem_dim(dest).x <= dest_p.x  ) return;

        T* dd = ui_c_updated(dest)(0,0);
        T* sd = ui_c_current(src)(0,0);

        for(size_t x = 0; x < size.x; x++)
            PTF(dest,dd,dim2(dest_p.x + x, dest_p.y),
                src, sd,dim2(src_p.x + x,  src_p.y),
                size.y, alfa);            
        __A_TIME_STOP
    }

    template<typename T, void(*PTF)(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 dpos,
                                    const maquis::types::p_dense_matrix_impl<T>& src, T *sd, dim2 spos, 
                                    size_t w, T alfa)>
    inline void __a_memptf(maquis::types::p_dense_matrix_impl<T>& dest, dim2 dest_p, 
                           const maquis::types::p_dense_matrix_impl<T>& src, dim2 src_p, 
                           dim2 size, T alfa = 0.0)
    {
        __A_TIME("ambient_memptf_f_kernel");
        // the ouput (dest) must be a pinned p_dense_matrix

        T* dd = ui_c_updated(dest)(ctxt.get_block_id().x, ctxt.get_block_id().y);
        size_t dx = ctxt.get_block_id().x * ui_c_get_mem_dim(dest).x;
        size_t dy = ctxt.get_block_id().y * ui_c_get_mem_dim(dest).y;
    
#ifdef AMBIENT_CHECK_BOUNDARIES
        if(ui_c_get_dim(dest).x - dest_p.x < size.x || ui_c_get_dim(dest).y - dest_p.y < size.y ||
           ui_c_get_dim(src).x - src_p.x   < size.x || ui_c_get_dim(src).y - src_p.y   < size.y) 
            maquis::cout << "Error: invalid memory movement" << std::endl;
#endif
    
        if( size.x == 0 || size.y == 0            ||
            dy + ui_c_get_mem_dim(dest).y <= dest_p.y  || 
            dx + ui_c_get_mem_dim(dest).x <= dest_p.x  ||
            dy >= dest_p.y + size.y               || 
            dx >= dest_p.x + size.x               ) return;

    // lets find dest-block starting point
        size_t offsetx = (dx < dest_p.x) ? dest_p.x % ui_c_get_mem_dim(dest).x : 0; dx += offsetx;
        size_t offsety = (dy < dest_p.y) ? dest_p.y % ui_c_get_mem_dim(dest).y : 0; dy += offsety;

        size_t sx = dx - dest_p.x + src_p.x;
        size_t sy = dy - dest_p.y + src_p.y;
    // lets find dest-block copy limits
        size_t sizey = __a_get_limit_y(dest, (dest_p.y + size.y)) - offsety;
        size_t sizex = __a_get_limit_x(dest, (dest_p.x + size.x)) - offsetx;

    // how many blocks do we need for this one
        int remaining = (int)sizey - (int)(ui_c_get_mem_dim(src).y - sy % ui_c_get_mem_dim(src).y);
        size_t num_src_blocks = remaining > 0 ? __a_ceil( remaining / ui_c_get_mem_dim(src).y ) + 1 : 1;

        dim2 dpos,spos; size_t v, w;
        for(size_t x = 0; x < sizex; x++){
            w = sizey;
            dpos.y = offsety;
            dpos.x = offsetx + x;
            spos.y = sy % ui_c_get_mem_dim(src).y;
            spos.x = (sx+x) % ui_c_get_mem_dim(src).x;
            for(int k = 0; k < num_src_blocks; k++){
                T* sd = ui_c_current(src)((sx+x) / ui_c_get_mem_dim(src).x, sy / ui_c_get_mem_dim(src).y + k);
                v = std::min(w, ui_c_get_mem_dim(src).y-spos.y);
                PTF(dest,dd,dpos,src,sd,spos,v,alfa);            
                w -= v;
                dpos.y += v;
                spos.y = 0;
            }
        }
        __A_TIME_STOP
    }

    template<typename T, void(*PTF)(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, dim2 dpos,
                                    const maquis::types::p_dense_matrix_impl<T>& src, T *sd, dim2 spos, 
                                    size_t w, T alfa)>
    inline void __a_memptf_reverse(maquis::types::p_dense_matrix_impl<T>& dest, dim2 dest_p, 
                                   const maquis::types::p_dense_matrix_impl<T>& src, dim2 src_p, 
                                   dim2 size, T alfa = 0.0)
    {
        __A_TIME("ambient_memptf_revers_f_kernel");
        // the input (src) must be a pinned p_dense_matrix

        dim2 context(ctxt.get_block_id().x, ctxt.get_block_id().y);

        T* sd = ui_c_current(src)(context.x, context.y);
        size_t sx = context.x * ui_c_get_mem_dim(src).x;
        size_t sy = context.y * ui_c_get_mem_dim(src).y;
    
#ifdef AMBIENT_CHECK_BOUNDARIES
        if(ui_c_get_dim_u(dest).x - dest_p.x < size.x || ui_c_get_dim_u(dest).y - dest_p.y < size.y ||
           ui_c_get_dim(src).x - src_p.x   < size.x || ui_c_get_dim(src).y - src_p.y   < size.y) 
            maquis::cout << "Error: invalid memory movement" << std::endl;
#endif
    
        if( size.x == 0 || size.y == 0          ||
            sy + ui_c_get_mem_dim(src).y <= src_p.y  || 
            sx + ui_c_get_mem_dim(src).x <= src_p.x  ||
            sy >= src_p.y + size.y              || 
            sx >= src_p.x + size.x               ) return;

    // lets find src-block starting point
        size_t offsety = (sy < src_p.y) ? src_p.y % ui_c_get_mem_dim(src).y : 0; sy += offsety;
        size_t offsetx = (sx < src_p.x) ? src_p.x % ui_c_get_mem_dim(src).x : 0; sx += offsetx;

        size_t dy = sy - src_p.y + dest_p.y;
        size_t dx = sx - src_p.x + dest_p.x;
    // lets find src-block copy limits
        size_t sizey = __a_get_limit_y(src, (src_p.y + size.y)) - offsety;
        size_t sizex = __a_get_limit_x(src, (src_p.x + size.x)) - offsetx;

    // how many blocks do we need for this one
        int remaining = (int)sizey - (int)(ui_c_get_mem_dim_u(dest).y - dy % ui_c_get_mem_dim_u(dest).y);
        size_t num_dest_blocks = remaining > 0 ? __a_ceil( remaining / ui_c_get_mem_dim_u(dest).y ) + 1 : 1;

        dim2 dpos,spos; size_t v, w;
        for(size_t x = 0; x < sizex; x++){
            w = sizey;
            spos.y = offsety;
            spos.x = offsetx + x;
            dpos.y = dy % ui_c_get_mem_dim_u(dest).y;
            dpos.x = (dx+x) % ui_c_get_mem_dim_u(dest).x;
            for(int k = 0; k < num_dest_blocks; k++){
                T* dd = ui_c_updated(dest)((dx+x) / ui_c_get_mem_dim_u(dest).x, dy / ui_c_get_mem_dim_u(dest).y + k);
                v = std::min(w, ui_c_get_mem_dim_u(dest).y-dpos.y);
                PTF(dest,dd,dpos,src,sd,spos,v,alfa);            
                w -= v;
                spos.y += v;
                dpos.y = 0;
            }
        }
        __A_TIME_STOP
    }

    /*template<void(*ASSIGN)(const models::v_model::object&, int, int), typename T>
    inline void block_2d_cycle(T& target){
    ///////////////////////////////////////////// 2D-block-cyclic decomposition
        int np = ctxt.np = 1; // can be a function arg   // process grid's num of rows 
        int nq = ctxt.nq = (int)(ctxt.get_size() / np); // process grid's num of cols 
        int rank_i = (int)(ctxt.get_rank() / nq); // process row
        int rank_j = (int)(ctxt.get_rank() % nq); // process col
    ///////////////////////////////////////////////////////////////////////////
        size_t sizey = get_grid_dim(target).y;
        size_t sizex = get_grid_dim(target).x;
        for(int i = rank_i; i < sizey; i += np){
            for(int j = rank_j; j < sizex; j += nq){
                ASSIGN(target, i, j);
            }
        }
    }*/

    template <typename T>
    inline void __a_atomic_refresh(maquis::types::p_dense_matrix_impl<T>& m){
        T* dm = ui_c_current(m)(0,0);
        T* rm = ui_c_updated(m)(0,0);
        __a_copy(rm, dm, ui_c_get_mem_dim(m).square());
    }

    inline void atomic_assign(revision& r, int x, int y){ 
        assign(r, x, y); 
    }

    inline void atomic_pin(revision& r, int x, int y){ 
        cfunctor* o = ctxt.get_op();
        o->add_condition();
        pin(o, r, x, y);
        assign(r, x, y);
    }

    inline void block_outright_assign(revision& r){
        size_t sizex = ui_l_get_grid_dim(r).x;
        size_t sizey = ui_l_get_grid_dim(r).y;
        for(int x = 0; x < sizex; x++)
        for(int y = 0; y < sizey; y++)
        assign(r, x, y);
    }

    inline void block_outright_pin(revision& r){
        cfunctor* o = ctxt.get_op();
        size_t sizex = ui_l_get_grid_dim(r).x;
        size_t sizey = ui_l_get_grid_dim(r).y;
        o->add_condition(sizey*sizex);
        for(int x = 0; x < sizex; x++)
        for(int y = 0; y < sizey; y++)
        pin(o, r, x, y);

        block_outright_assign(r);
    }

    inline void block_2d_cycle_assign(revision& r){ block_outright_assign(r); }
    inline void block_2d_cycle_pin(revision& r){ block_outright_pin(r); }
    inline void block_outright_conditional_assign(revision& r){ block_outright_pin(r); }
    inline void block_2d_cycle_conditional_assign(revision& r){ block_outright_conditional_assign(r); }

}

#endif
