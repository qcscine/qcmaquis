#ifndef VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP

#include "kernels_cpu.hpp"

namespace vli {

    /**
    template forward declaration 
    */
    template <class BaseInt, int Size>
    class vli_cpu;

namespace detail {


template <class BaseInt, int Size>
void plus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b )
{
    addition_classic_cpu<BaseInt, Size>(&vli_a[0],&vli_b[0]);
}


template <class BaseInt, int Size>
void multiplies_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b )
{
     multiplication_classic_cpu<BaseInt,Size>(&vli_a[0], &vli_b[0]);
}

} //namespace detail
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
