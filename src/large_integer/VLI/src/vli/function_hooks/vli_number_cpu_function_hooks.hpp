#ifndef VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP


//#include "vli/detail/kernels_cpu_gpu.hpp"
#include "vli/detail/kernels_cpu_asm.h"

namespace vli 
{

    /**
    template forward declaration 
    */
    template <class BaseInt, std::size_t Size>
    class vli_cpu;

    template <class BaseInt, std::size_t Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b ){
        detail::add192(&vli_a[0],&vli_b[0]);
    }
     
    template <class BaseInt, std::size_t Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b ){  
        assert(b < max_value<BaseInt>::value); //avoid overflow
        detail::add64(&vli_a[0],&b);
    }
    
    template <class BaseInt, std::size_t Size>
    void minus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b ){
        detail::sub192(&vli_a[0],&vli_b[0]);
    }
    
    template <class BaseInt, std::size_t Size>
    void minus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b ){  
        assert(b < max_value<BaseInt>::value); //avoid overflow
        detail::sub64(&vli_a[0],&b);
    }

    template <class BaseInt, std::size_t Size>
    void multiplies(vli_cpu<BaseInt, 2*Size>& vli_res , vli_cpu<BaseInt,Size> const & vli_a, vli_cpu<BaseInt,Size> const & vli_b){

        //      detail::kernels_multiplication_classic<BaseInt,Size>(&vli_res[0],&vli_a[0], &vli_b[0]);
    }
    
    template <class BaseInt, std::size_t Size>
    void multiplies_assign( vli_cpu<BaseInt, Size>& vli_a , vli_cpu<BaseInt,Size> const & vli_b){ 
        detail::mul192(&vli_a[0],&vli_b[0]);
    }

    template <class BaseInt, std::size_t Size>
    void multiplies_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b){
        detail::mul64(&vli_a[0],&b);
    }
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
