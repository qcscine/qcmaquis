#include "kernels_gpu.h"

#ifndef VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP

namespace vli {

    /**
    template forward declaration 
    */
    template <class BaseInt>
    class vli_gpu;

namespace detail {

template <class BaseInt>
void plus_assign(vli_gpu<BaseInt> & vli_a, vli_gpu<BaseInt> const& vli_b )
{
    plus_assign_gpu( vli_a.p(), vli_b.p(), 1, vli_gpu<BaseInt>::size);
}


template <class BaseInt>
void multiplies_assign(vli_gpu<BaseInt> & vli_a, vli_gpu<BaseInt> const& vli_b )
{
    // TODO this seems to be quite inefficient
    // (We implement the a *= b as a regular c = a * b and do a swap(a,c)...)
    // (I guess it's better to try implement a *= b directly.)
    vli_gpu<BaseInt> vli_c;
    entrywise_multiplies_gpu( vli_a.p(), vli_b.p(), vli_c.p(), 1, vli_gpu<BaseInt>::size);
    swap(vli_a,vli_c);
}

} //namespace detail
} //namespace vli

#endif //VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
