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

    template<typename V, typename T>
    inline void* __a_solidify_atomic(const iteratable<T>& o){
        c_revision& r = o.ui_c_revision_0();
        void* memory = malloc(o.spec.size);
        memcpy(memory, (char*)r, o.spec.size);
        return memory;
    }

    template<typename V, typename T>
    inline void __a_disperse_atomic(void* data, iteratable<T>& o){
        w_revision& r = o.ui_w_revision_1();
        memcpy((char*)r, data, o.spec.size);
        free(data);
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
        __a_copy(&dd[dpos.x*ui_c_get_dim(dest).y+dpos.y],
                 &sd[spos.x*ui_c_get_dim(src).y+spos.y],
                 w);
    }

    template <typename T>
    inline void __a_memcpy(T* dd, T *sd, size_t w, T alfa){
        memcpy(dd, sd, w);
    }

    template <typename T>
    inline void __a_memscal(matrix_impl<T>& dest, T* dd, dim2 dpos, const matrix_impl<T>& src, T *sd, dim2 spos, size_t w, T alfa){
        for(int z = 0; z < w; z++)
            dd[dpos.x*ui_c_get_dim(dest).y+dpos.y + z] += sd[spos.x*ui_c_get_dim(src).y+spos.y + z]*alfa; // be carefull that dd != sd
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
        int lda = ui_c_get_dim(src).y;
        int ldb = ui_c_get_dim(dst).y;

        T* sd = (T*)ui_c_current(src) + src_p.y + src_p.x*lda;
        T* dd = (T*)ui_r_updated(dst) + dst_p.y + dst_p.x*ldb;

        do{ PTF(dd, sd, m, alfa); sd += lda; dd += ldb; }while(--n > 0);
        __A_TIME_C_STOP
    }

    template <typename T>
    inline void __a_atomic_refresh(matrix_impl<T>& m){
        T* dm = ui_c_current(m);
        T* rm = ui_w_updated(m);
        if(dm != rm) memcpy(rm, dm, ui_c_get_mem_size(m));
    }

} } }

#endif
