#ifndef __MAQUIS_TYPES_KERNELS_UTILS_HPP__
#define __MAQUIS_TYPES_KERNELS_UTILS_HPP__

#include <limits>
#include "utils/timings.h"

//#define AMBIENT_COMPUTATIONAL_TIMINGS
//#define AMBIENT_CHECK_BOUNDARIES

#ifdef AMBIENT_CHECK_BOUNDARIES
#include <execinfo.h>
#endif

extern "C" {
    double ddot_(const int*, const double*, const int*, const double*, const int*);
}

namespace ambient { namespace numeric { namespace kernels {

    using ambient::numeric::matrix_impl;

    #include "ambient/utils/numeric.h" // BLAS/LAPACK prototypes
    #include "ambient/utils/ceil.h"
   
    #ifdef AMBIENT_COMPUTATIONAL_TIMINGS
        #define __A_TIME_C(name) static __a_timer time(name); time.begin();
        #define __A_TIME_C_STOP time.end();
    #else
        #define __A_TIME_C(name) 
        #define __A_TIME_C_STOP 
    #endif

    // {{{ continuous memory mangling
    template<typename V, typename T>
    inline void* __a_solidify(const iteratable<T>& o){
        using ambient::models::velvet::memspec;

        c_revision& r = o.ui_c_revision_0();
        size_t iterator = 0;
        char* memory = NULL;
        size_t stride = o.spec.block.y*sizeof(V);
        size_t block = o.spec.block.x;
        dim2 grid = o.spec.grid;

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
        using ambient::models::velvet::memspec;
        using ambient::models::velvet::history;

        w_revision& r = o.ui_w_revision_1();
        char* memory = (char*)data;
        size_t stride = o.spec.block.y*sizeof(V);
        size_t block = o.spec.block.x;
        dim2 grid = o.spec.grid;

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

    template<typename V, typename T>
    inline void* __a_solidify_atomic(const iteratable<T>& o){
        c_revision& r = o.ui_c_revision_0();
        void* memory = malloc(o.spec.size);
        memcpy(memory, (char*)r(0,0), o.spec.size);
        return memory;
    }

    template<typename V, typename T>
    inline void __a_disperse_atomic(void* data, iteratable<T>& o){
        w_revision& r = o.ui_w_revision_1();
        memcpy((char*)r(0,0), data, o.spec.size);
        free(data);
    }
    // }}}

    template<typename T>
    inline size_t __a_get_limit_x(const T& a, size_t n){
        return std::min(ui_c_get_mem_dim(a).x, n-ctxt.get_block_id().x*ui_c_get_mem_dim(a).x);
    }

    template<typename T>
    inline size_t __a_get_limit_x(const T& a){
        return std::min(ui_c_get_mem_dim(a).x, ui_c_get_dim(a).x-ctxt.get_block_id().x*ui_c_get_mem_dim(a).x);
    }

    template<typename T>
    inline size_t __a_get_limit_y(const T& a, size_t m){
        return std::min(ui_c_get_mem_dim(a).y, m-ctxt.get_block_id().y*ui_c_get_mem_dim(a).y);
    }

    template<typename T>
    inline size_t __a_get_limit_y(const T& a){
        return std::min(ui_c_get_mem_dim(a).y, ui_c_get_dim(a).y-ctxt.get_block_id().y*ui_c_get_mem_dim(a).y);
    }

    template <typename T> inline T __a_dot(T* a, T* b, int size){
        T summ(0);
        for(size_t k=0; k < size; k++)
           summ += a[k]*b[k];
        return summ;
    }

    inline double __a_dot(double* a, double* b, int size){
        static const int ONE = 1;
        return ddot_(&size, a, &ONE, b, &ONE);
    }

    template<typename T>
    inline void __a_copy(T* dst, T* src, int size){
        memcpy(dst, src, size*sizeof(T));
    }

    template <typename T>
    inline void __a_memcpy(matrix_impl<T>& dest, T* dd, dim2 dpos, const matrix_impl<T>& src, T *sd, dim2 spos, size_t w, T alfa){
        __a_copy(&dd[dpos.x*ui_c_get_mem_dim(dest).y+dpos.y],
                 &sd[spos.x*ui_c_get_mem_dim(src).y+spos.y],
                 w);
    }

    template <typename T>
    inline void __a_memcpy(T* dd, T *sd, size_t w, T alfa){
        memcpy(dd, sd, w);
    }

    template <typename T>
    inline void __a_memscal(matrix_impl<T>& dest, T* dd, dim2 dpos, const matrix_impl<T>& src, T *sd, dim2 spos, size_t w, T alfa){
        for(int z = 0; z < w; z++)
            dd[dpos.x*ui_c_get_mem_dim(dest).y+dpos.y + z] += sd[spos.x*ui_c_get_mem_dim(src).y+spos.y + z]*alfa; // be carefull that dd != sd
    }

    template <typename T>
    inline void __a_memscal(T* dd, T *sd, size_t w, T alfa){
        int z = w/sizeof(T);
        do{ *dd++ += alfa*(*sd++); }while(--z > 0); // be carefull that dd != sd
    }

    template<typename T, void(*PTF)(T* dd, T* sd, size_t w, T alfa)>
    inline void __a_memptf_atomic_r(matrix_impl<T>& dst, dim2 dst_p, 
                                    const matrix_impl<T>& src, dim2 src_p, 
                                    dim2 size, T alfa = 0.0)
    {
        __A_TIME_C("ambient_memptf_fr_atomic_kernel");
#ifdef AMBIENT_CHECK_BOUNDARIES
        if(ui_c_get_dim(dst).x - dst_p.x < size.x || ui_c_get_dim(dst).y - dst_p.y < size.y ||
           ui_c_get_dim(src).x - src_p.x < size.x || ui_c_get_dim(src).y - src_p.y < size.y){
            ambient::cout << "Error: invalid memory movement: " << std::endl;
            ambient::cout << "Matrix dst " << ui_c_get_dim(dst).x << "x" << ui_c_get_dim(dst).y << "\n";
            ambient::cout << "Dest p " << dst_p.x << "x" << dst_p.y << "\n";
            ambient::cout << "Matrix src " << ui_c_get_dim(src).x << "x" << ui_c_get_dim(src).y << "\n";
            ambient::cout << "Src p " << src_p.x << "x" << src_p.y << "\n";
            ambient::cout << "Block size " << size.x << "x" << size.y << "\n";

            void *array[10];
            size_t size = backtrace(array, 10);
            backtrace_symbols_fd(array, size, 2);
        }
#endif
        int n = size.x;
        int m = size.y*sizeof(T);
        int lda = ui_c_get_mem_dim(src).y;
        int ldb = ui_c_get_mem_dim(dst).y;

        T* sd = (T*)ui_c_current(src)(0,0) + src_p.y + src_p.x*lda;
        T* dd = (T*)ui_r_updated(dst)(0,0) + dst_p.y + dst_p.x*ldb;

        do{ PTF(dd, sd, m, alfa); sd += lda; dd += ldb; }while(--n > 0);
        __A_TIME_C_STOP
    }

    template<typename T, void(*PTF)(matrix_impl<T>& dest, T* dd, dim2 dpos,
                                    const matrix_impl<T>& src, T *sd, dim2 spos, 
                                    size_t w, T alfa)>
    inline void __a_memptf(matrix_impl<T>& dest, dim2 dest_p, 
                           const matrix_impl<T>& src, dim2 src_p, 
                           dim2 size, T alfa = 0.0)
    {
        __A_TIME_C("ambient_memptf_f_kernel");
        // the ouput (dest) must be a pinned matrix

        T* dd = ui_r_updated(dest)(ctxt.get_block_id().x, ctxt.get_block_id().y);
        size_t dx = ctxt.get_block_id().x * ui_c_get_mem_dim(dest).x;
        size_t dy = ctxt.get_block_id().y * ui_c_get_mem_dim(dest).y;
    
#ifdef AMBIENT_CHECK_BOUNDARIES
        if(ui_c_get_dim(dest).x - dest_p.x < size.x || ui_c_get_dim(dest).y - dest_p.y < size.y ||
           ui_c_get_dim(src).x - src_p.x   < size.x || ui_c_get_dim(src).y - src_p.y   < size.y) 
            ambient::cout << "Error: invalid memory movement" << std::endl;
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
        __A_TIME_C_STOP
    }

    template<typename T, void(*PTF)(matrix_impl<T>& dest, T* dd, dim2 dpos,
                                    const matrix_impl<T>& src, T *sd, dim2 spos, 
                                    size_t w, T alfa)>
    inline void __a_memptf_reverse(matrix_impl<T>& dest, dim2 dest_p, 
                                   const matrix_impl<T>& src, dim2 src_p, 
                                   dim2 size, T alfa = 0.0)
    {
        __A_TIME_C("ambient_memptf_revers_f_kernel");
        // the input (src) must be a pinned matrix

        dim2 context(ctxt.get_block_id().x, ctxt.get_block_id().y);

        T* sd = ui_c_current(src)(context.x, context.y);
        size_t sx = context.x * ui_c_get_mem_dim(src).x;
        size_t sy = context.y * ui_c_get_mem_dim(src).y;
    
#ifdef AMBIENT_CHECK_BOUNDARIES
        if(ui_c_get_dim(dest).x - dest_p.x < size.x || ui_c_get_dim(dest).y - dest_p.y < size.y ||
           ui_c_get_dim(src).x - src_p.x   < size.x || ui_c_get_dim(src).y - src_p.y   < size.y) 
            ambient::cout << "Error: invalid memory movement" << std::endl;
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
        int remaining = (int)sizey - (int)(ui_c_get_mem_dim(dest).y - dy % ui_c_get_mem_dim(dest).y);
        size_t num_dest_blocks = remaining > 0 ? __a_ceil( remaining / ui_c_get_mem_dim(dest).y ) + 1 : 1;

        dim2 dpos,spos; size_t v, w;
        for(size_t x = 0; x < sizex; x++){
            w = sizey;
            spos.y = offsety;
            spos.x = offsetx + x;
            dpos.y = dy % ui_c_get_mem_dim(dest).y;
            dpos.x = (dx+x) % ui_c_get_mem_dim(dest).x;
            for(int k = 0; k < num_dest_blocks; k++){
                T* dd = ui_r_updated(dest)((dx+x) / ui_c_get_mem_dim(dest).x, dy / ui_c_get_mem_dim(dest).y + k);
                v = std::min(w, ui_c_get_mem_dim(dest).y-dpos.y);
                PTF(dest,dd,dpos,src,sd,spos,v,alfa);            
                w -= v;
                spos.y += v;
                dpos.y = 0;
            }
        }
        __A_TIME_C_STOP
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
    inline void __a_atomic_refresh(matrix_impl<T>& m){
        T* dm = ui_c_current(m)(0,0);
        T* rm = ui_w_updated(m)(0,0);
        if(dm != rm) memcpy(rm, dm, ui_c_get_mem_size(m));
    }

} } }

#endif
