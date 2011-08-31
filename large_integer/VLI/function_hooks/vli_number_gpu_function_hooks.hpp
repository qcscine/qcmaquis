
#ifndef VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP

#include "detail/kernels_gpu.h"


namespace vli {

/**
template forward declaration 
*/
template <class BaseInt, int Size>
class vli_gpu;


template <class BaseInt, int Size>
void plus_assign(vli_gpu<BaseInt,Size> & vli_a, vli_gpu<BaseInt,Size> const& vli_b)
{
    detail::plus_assign(vli_a.p(), vli_b.p());
}

template <class BaseInt, int Size>
void plus_assign(vli_gpu<BaseInt,Size> & vli_a, BaseInt b)
{
//    assert(b < max_value<BaseInt>::value); //avoid overflow
    detail::plus_assign(vli_a.p(), b);
}

template <class BaseInt, int Size>
void minus_assign(vli_gpu<BaseInt,Size> & vli_a, vli_gpu<BaseInt,Size> const& vli_b )
{
    detail::minus_assign(vli_a.p(),vli_b.p());
}

template <class BaseInt, int Size>
void multiplies_assign(vli_gpu<BaseInt,Size>& vli_a, BaseInt a)
{
    detail::multiplies_assign(vli_a.p(),a);
}

template <class BaseInt, int Size>
void multiplies_assign(vli_gpu<BaseInt,Size> const& vli_a, vli_gpu<BaseInt,Size> const& vli_b, vli_gpu<BaseInt, Size> & vli_c)
{
    detail::entrywise_multiplies( vli_a.p(), vli_b.p(), vli_c.p());
}

} //namespace vli

#endif //VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
