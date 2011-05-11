#ifndef VLI_VECTOR_GPU_FUNCTION_HOOKS_HPP
#define  VLI_VECTOR_GPU_FUNCTION_HOOKS_HPP

#include "kernels_gpu.h"

namespace vli {

    template <class BaseInt>
    class vli_vector_gpu;

namespace detail {

template <class BaseInt>
void plus_assign(vli_vector_gpu<BaseInt>& v_a, vli_vector_gpu<BaseInt> const& v_b)
{
    assert( v_a.size() == v_b.size() );
    plus_assign_gpu(v_a.p(), v_b.p(), v_a.size(), vli_vector_gpu<BaseInt>::vli_size);
}

template <class BaseInt>
void multiplies_assign(vli_vector_gpu<BaseInt>& v_a, vli_gpu<BaseInt> const& v)
{
    // NOT IMPLEMENTED YET
    assert(false);
}

template <class BaseInt>
void entrywise_multiplies_assign(vli_vector_gpu<BaseInt>& v_a, vli_vector_gpu<BaseInt> const& v_b)
{
    assert( v_a.size() == v_b.size() );
    entrywise_multiplies_assign_gpu(v_a.p(), v_b.p(), v_a.size(), vli_vector_gpu<BaseInt>::vli_size);
}

} //namespace detail
} //namespace vli

#endif // VLI_VECTOR_GPU_FUNCTION_HOOKS_HPP
