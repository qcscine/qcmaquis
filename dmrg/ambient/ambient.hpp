#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include "ambient/ambient.h"

namespace blas{ 
//  forward declarations if handy
    template <typename T, ambient::policy P = ambient::ANY> 
    class p_dense_matrix; 
}

namespace ambient{ 

// synchronization primitives //
    void memoryfence();

// kernels management //
    template <typename FL, typename FC, class T0, class T1, class T2>
    void push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1, T2& arg2);

    template <typename ST, typename FL, typename FC, class T0, class T1>
    ST& push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1);

    template <typename L, typename R>
    void pin(L& proxy_object, const R& real_object);

    template<typename T>
    void assign(T& ref, int i, int j = 0);

    template<typename T>
    void assign(const T& ref, int i, int j = 0);

// breakdown information //
    template<typename T>
    inline std::pair<unsigned int*,size_t> get_id(T& ref);

    template<typename T>
    inline dim2 get_dim(T& ref);

    template<typename T>
    inline dim2 get_work_dim(T& ref);

    template<typename T>
    inline dim2 get_gpu_dim(T& ref);

    template<typename T>
    inline dim2 get_grid_dim(T& ref);

    template<typename T>
    inline dim2 get_mem_dim(T& ref);

    template<typename T>
    inline dim2 get_mem_t_dim(T& ref);

    template<typename T>
    inline dim2 get_item_dim(T& ref);

    template<typename T>
    inline dim2 get_block_id(T& ref);

    #include "ambient/interface/core.hpp"
    #include "ambient/interface/i_kernels.hpp"
    #include "ambient/interface/l_kernels.hpp"
    #include "ambient/interface/c_kernels.hpp"
    #include "ambient/interface/c_scalapack_kernels.hpp"
    #include "ambient/interface/model.hpp"
    #include "ambient/interface/interface.hpp"
}
#endif
