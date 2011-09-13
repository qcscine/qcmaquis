#ifndef VLI_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP
#define VLI_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP

#include "vli/detail/kernels_gpu.h"

namespace vli{
    
    using detail::vli_size_tag;

    template<class Vli, int Order>
    class polynomial_gpu;

    template <class BaseInt, int Size, int Order>
    void poly_addition_int(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & p1, int b)
    {
        detail::plus_assign_poly_int(vli_size_tag<Size>(),p1.p(),b);
    } 

    template <class BaseInt, int Size, int Order, class Monomial> // monomial can be a monomial_gpu or a vli_gpu
    void poly_mono_multiply(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                            polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p1, 
                            Monomial const & m2)
    {
        detail::poly_mono_multiply(vli_size_tag<Size>(),p1.p(),m2.p(),result.p(), m2.j_exp_, m2.h_exp_);
    }
    
    template <class BaseInt, int Size, int Order>
    void plus_assign_poly(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p)
    {
        detail::plus_assign_poly(vli_size_tag<Size>(),result.p(),p.p());
    }

    template <class BaseInt, int Size, int Order>
    void poly_substraction(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p)
    {
        detail::minus_assign_poly(vli_size_tag<Size>(),result.p(),p.p());
    }

    template <class BaseInt, int Size, int Order>
    void poly_multiply(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p1, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p2)
    {
        detail::poly_poly_multiply(vli_size_tag<Size>(),p1.p(),p2.p(),result.p());
    }
}// end namespace 

#endif
