#ifndef VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP


#include "vli/detail/kernels_cpu_gpu.hpp"

namespace vli 
{

    /**
    template forward declaration 
    */
    template <class BaseInt, int Size>
    class vli_cpu;

    template <class BaseInt, int Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b )
    {
        detail::kernels_addition_classic<BaseInt, Size>(&vli_a[0],&vli_b[0]);    
    }
     
    template <class BaseInt, int Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const& b )
    {  
        assert(b < max_value<BaseInt>::value); //avoid overflow
        detail::kernels_addition_block<BaseInt, Size>(&vli_a[0],&b);  
    }
    
    template <class BaseInt, int Size>
    void multiplies_assign( vli_cpu<BaseInt, Size>& vli_a , vli_cpu<BaseInt,Size> const & vli_b)
    { 
        vli_cpu<BaseInt, Size> vli_res;
        detail::kernels_multiplication_classic_truncate<BaseInt,Size>(&vli_res[0],&vli_a[0], &vli_b[0]);
        vli_a = vli_res;
    }
    
    template <class BaseInt, int Size>
    void multiplies_assign(vli_cpu<BaseInt,Size> & vli_a, int n)
    {
        BaseInt tmp = static_cast<BaseInt>(n);
        detail::kernels_multiplication_number<BaseInt,Size>(&vli_a[0],tmp); 
    }
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
