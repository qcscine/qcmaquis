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
    entrywise_multiplies_assign_gpu( vli_a.p(), vli_b.p(), 1, vli_gpu<BaseInt>::size);
}

} //namespace detail
} //namespace vli

#endif //VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
