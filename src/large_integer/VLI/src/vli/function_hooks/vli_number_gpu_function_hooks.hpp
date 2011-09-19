
#ifndef VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP

#include "vli/detail/kernels_gpu.h"


namespace vli {

/**
template forward declaration 
*/
template <class BaseInt, std::size_t Size>
class vli_gpu;


template <class BaseInt, std::size_t Size>
void plus_assign(vli_gpu<BaseInt,Size> & vli_a, vli_gpu<BaseInt,Size> const& vli_b)
{
    using detail::vli_size_tag;
    detail::plus_assign(vli_size_tag<Size>(), vli_a.p(), vli_b.p());
}

template <class BaseInt, std::size_t Size>
void plus_assign(vli_gpu<BaseInt,Size> & vli_a, BaseInt b)
{
//    assert(b < max_value<BaseInt>::value); //avoid overflow
    using detail::vli_size_tag;
    detail::plus_assign(vli_size_tag<Size>(), vli_a.p(), b);
}

template <class BaseInt, std::size_t Size>
void minus_assign(vli_gpu<BaseInt,Size> & vli_a, vli_gpu<BaseInt,Size> const& vli_b )
{
    using detail::vli_size_tag;
    vli_gpu<BaseInt,Size> tmp(vli_b);
    detail::minus_assign_destructive(vli_size_tag<Size>(), vli_a.p(),tmp.p());
}

template <class BaseInt, std::size_t Size>
void minus_assign(vli_gpu<BaseInt,Size> & vli_a, int a )
{
    using detail::vli_size_tag;
    vli_gpu<BaseInt,Size> tmp(a);
    detail::minus_assign_destructive(vli_size_tag<Size>(), vli_a.p(),tmp.p());
}

template <class BaseInt, std::size_t Size>
void multiplies_assign(vli_gpu<BaseInt,Size>& vli_a, BaseInt a)
{
    using detail::vli_size_tag;
    detail::multiplies_assign(vli_size_tag<Size>(), vli_a.p(), a);
}

template <class BaseInt, std::size_t Size>
void multiplies_assign(vli_gpu<BaseInt,Size> const& vli_a, vli_gpu<BaseInt,Size> const& vli_b, vli_gpu<BaseInt, Size> & vli_c)
{
    using detail::vli_size_tag;
    detail::entrywise_multiplies(vli_size_tag<Size>(), vli_a.p(), vli_b.p(), vli_c.p());
}

} //namespace vli

#endif //VLI_NUMBER_GPU_FUNCTION_HOOKS_HPP
