#ifndef VLI_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP
#define VLI_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP

#include "kernels_gpu.h"

namespace vli{
    
    template<class Vli, int Order>
    class polynomial_gpu;
    /*
    template <class BaseInt, int Size, int Order>
    void poly_multiply(polynomial<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial<vli_gpu<BaseInt, Size>, Order> const & p1, 
                       polynomial<vli_gpu<BaseInt, Size>, Order> const & p2)
    {
        typedef typename polynomial<vli_gpu<BaseInt, Size>, Order>::size_type size_type;	
        enum {vli_size  = polynomial<vli_gpu<BaseInt, Size>, Order>::vli_size}; 
        enum {max_order = polynomial<vli_gpu<BaseInt, Size>, Order>::max_order}; 
        detail::poly_multiply_gpu(p1(0,0).p(),p2(0,0).p(),result(0,0).p(),vli_size,max_order);
    }
    */
    template <class BaseInt, int Size, int Order>
    void poly_multiply(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p1, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p2)
    {
        //C
        //C In principle enums were correct, but why don't we just use the template parameters directly?
        //C
        detail::poly_multiply_gpu(p1.p(),p2.p(),result.p(),Size,Order);
    }
    
    template <class BaseInt, int Size, int Order>
    void poly_addition(polynomial_gpu<vli_gpu<BaseInt, Size>, Order> & result, 
                       polynomial_gpu<vli_gpu<BaseInt, Size>, Order> const & p)
    {
        //C
        //C In principle enums were correct, but why don't we just use the template parameters directly?
        //C
        detail::poly_addition_gpu(result.p(),p.p(),Size,Order);
    }
        
    

}// end namespace 

#endif
