#ifndef VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP

#include "kernels_cpu.hpp"

namespace vli {
    
    /**
     template forward declaration 
     */
    template <class BaseInt>
    class vli_cpu;
    
    namespace detail {
        
    template <class BaseInt>
    void multiplies_assign(vli_cpu<BaseInt> & vli_a, vli_cpu<BaseInt> const& vli_b )
    {
        multiplication_classic_cpu(&vli_a[0], &vli_b[0]);
    }
        
    } //namespace detail
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
