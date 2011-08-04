#ifndef VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP

#include "detail/kernels_cpu.hpp"

namespace vli {

    /**
    template forward declaration 
    */
    template <class BaseInt, int Size>
    class vli_cpu;

    template <class BaseInt, int Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b )
    {
        addition_classic_cpu<BaseInt, Size>(&vli_a[0],&vli_b[0]);
    }
     
    template <class BaseInt, int Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const& b )
    {  
        assert(b < max_value<BaseInt>::value);
        addition_classic_cpu<BaseInt, Size>(&vli_a[0],&b);
    }
    
    template <class BaseInt, int Size1, int Size2>
    void multiplies( vli_cpu<BaseInt, Size1> vli_res, vli_cpu<BaseInt,Size2> const & vli_a, vli_cpu<BaseInt,Size2> const& vli_b)
    {
        multiplication_classic_cpu<BaseInt,Size2>(&vli_res[0],&vli_a[0], &vli_b[0]);
    }
    
    template <class BaseInt, int Size>
    void multiplies_assign_array(BaseInt* a, BaseInt const* b )
    {
        multiplication_classic_cpu<BaseInt,Size>(a, b);
    }

    template <class BaseInt, int Size>
    void plus_assign_array(BaseInt* a, BaseInt const* b )
    {
        addition_classic_cpu<BaseInt,Size>(a, b);
    }
    
    template <class BaseInt, int Size>
    void multiplies_assign_number(vli_cpu<BaseInt,Size> & vli_a, int n)
    {
        BaseInt tmp = static_cast<BaseInt>(n);
        multiplication_classic_cpu_number<BaseInt,Size>(&vli_a[0],tmp);
    }
    
    
    
    
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
