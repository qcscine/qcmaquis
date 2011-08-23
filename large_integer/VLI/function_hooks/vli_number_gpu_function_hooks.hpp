#include "detail/kernels_gpu.h"

#ifndef VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP

namespace vli {

/**
template forward declaration 
*/
template <class BaseInt, int Size>
class vli_gpu;


template <class BaseInt, int Size>
void plus_assign(vli_gpu<BaseInt,Size> & vli_a, vli_gpu<BaseInt,Size> const& vli_b )
{
    using detail::plus_assign_gpu;
    plus_assign_lib(vli_a.p(), vli_b.p());
}

/*

template <class BaseInt, int Size>
void plus_assign(vli_gpu<BaseInt,Size> & vli_a, vli_gpu<BaseInt,Size> const& vli_b )
{
    using detail::plus_assign_gpu;
    plus_assign_gpu( vli_a.p(), vli_b.p(), Size);
}

*/
template <class BaseInt, int Size>
void plus_assign_int(vli_gpu<BaseInt,Size> & vli_a, BaseInt b )
{
    using detail::plus_assign_int_gpu;
    plus_assign_int_gpu( vli_a.p(), b, Size);
}

template <class BaseInt, int Size>
void minus_assign(vli_gpu<BaseInt,Size> & vli_a, vli_gpu<BaseInt,Size> const& vli_b )
{
    using detail::minus_assign_gpu;
    minus_assign_gpu(vli_a.p(),vli_b.p(),Size);
}

template <class BaseInt, int Size>
void multiplies_assign_single(vli_gpu<BaseInt,Size>& vli_a, BaseInt a)
{
    using detail::multiplies_assign_single_gpu;
    multiplies_assign_single_gpu(vli_a.p(),a,Size);

}

template <class BaseInt, int Size>
void multiplies_assign(vli_gpu<BaseInt,Size> const& vli_a, vli_gpu<BaseInt,Size> const& vli_b, vli_gpu<BaseInt, Size> & vli_c)
{
    using detail::entrywise_multiplies_gpu;
    entrywise_multiplies_gpu( vli_a.p(), vli_b.p(), vli_c.p(), Size);
}

} //namespace vli

#endif //VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
