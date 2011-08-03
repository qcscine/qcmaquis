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
    
    template <class BaseInt, int Size>
    void multiplies_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b )
    {
        multiplication_classic_cpu<BaseInt,Size>(&vli_a[0], &vli_b[0]);
    }

    template <class BaseInt, int Size>
    void multiplies_assign_negate(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b )
    {
        // no tmp variable, I minimize copy .....
        if (vli_a.is_negative()) // this is -
        {
            vli_a.negate(); // - to + 
            multiplication_classic_cpu<BaseInt,Size>(&vli_a[0], &vli_b[0]);
            vli_a.negate(); // + to -
        }
        else // vli_b is negative
        { // Putain de const -> const
            vli_cpu<BaseInt,Size>  tmp(vli_b); 
            tmp.negate();
            multiplication_classic_cpu<BaseInt,Size>(&vli_a[0], &tmp[0]);
            vli_a.negate();    
        }
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
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
