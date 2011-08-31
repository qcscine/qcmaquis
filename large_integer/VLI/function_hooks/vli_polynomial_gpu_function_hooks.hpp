#ifndef VLI_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP
#define VLI_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP

#include "detail/kernels_gpu.h"

namespace vli{
    
    template<class Vli, int Order>
    class polynomial_gpu;

    template <class BaseInt, int Size, int Order>
    void poly_addition_int(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & p1, int b)
    {
        detail::plus_assign_poly_int(p1.p(),b);
    } 

    template <class BaseInt, int Size, int Order, class Monomial> // monomial can be a monomial_gpu or a vli_gpu
    void poly_mono_multiply(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                            polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p1, 
                            Monomial const & m2)
    {
        detail::poly_mono_multiply(p1.p(),m2.p(),result.p());
    }
    
    template <class BaseInt, int Size, int Order>
    void plus_assign_poly(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p)
    {
        detail::plus_assign_poly(result.p(),p.p());
    }

    template <class BaseInt, int Size, int Order>
    void poly_substraction(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p)
    {
        detail::minus_assign_poly(result.p(),p.p());
    }

    template <class BaseInt, int Size, int Order>
    void poly_multiply(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p1, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p2)
    {
         detail::poly_poly_multiply(p1.p(),p2.p(),result.p());
    }
}// end namespace 

#endif
