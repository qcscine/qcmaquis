#ifndef VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP


//#include "vli/detail/kernels_cpu_gpu.hpp"
#include "vli/detail/kernels_cpu_asm.h"

namespace vli 
{
    //forwart declaration
    template <class BaseInt, std::size_t Size>
    class vli_cpu;
    
    //declaration wrapper
    //addition
    template<typename BaseInt, std::size_t Size>
    void add(BaseInt * x, BaseInt const b);
    
    template<typename BaseInt, std::size_t Size>
    void add(BaseInt * x, BaseInt const* y); 

    //substraction
    template<typename BaseInt, std::size_t Size>
    void sub(BaseInt * x, BaseInt const b);
    
    template<typename BaseInt, std::size_t Size>
    void sub(BaseInt * x, BaseInt const* y); 

    //multiplication
    template<typename BaseInt, std::size_t Size>
    void mul(BaseInt * x,BaseInt const y);

    template <class BaseInt, std::size_t Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b ){
        add<BaseInt,Size>(&vli_a[0],&vli_b[0]);
    }
     
    template <class BaseInt, std::size_t Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b ){  
        assert(b < max_value<BaseInt>::value); //avoid overflow
        add<BaseInt,Size>(&vli_a[0],b);
    }
    
    template <class BaseInt, std::size_t Size>
    void minus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b ){
        sub<BaseInt,Size>(&vli_a[0],&vli_b[0]);
    }
    
    template <class BaseInt, std::size_t Size>
    void minus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b ){  
        assert(b < max_value<BaseInt>::value); //avoid overflow
        sub<BaseInt,Size>(&vli_a[0],b);
    }

    template <class BaseInt, std::size_t Size>
    void multiplies(vli_cpu<BaseInt, 2*Size>& vli_res , vli_cpu<BaseInt,Size> const & vli_a, vli_cpu<BaseInt,Size> const & vli_b){
        detail::mul384_192_192(&vli_res[0],&vli_a[0],&vli_b[0]);
    }
    
    template <class BaseInt, std::size_t Size>
    void multiplies_assign( vli_cpu<BaseInt, Size>& vli_a , vli_cpu<BaseInt,Size> const & vli_b){ 
        detail::mul192_192(&vli_a[0],&vli_b[0]);
    }

    template <class BaseInt, std::size_t Size>
    void multiplies_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b){
        mul<BaseInt,Size>(&vli_a[0],b);
    }

    //specialization add
    template<>
    void add<unsigned long int,6>(unsigned long int * x,unsigned long int const* y){
        detail::add384_384(x,y); // 384 + 384 = 384
    }; 
    
    template<>
    void add<unsigned long int,3>(unsigned long int * x,unsigned long int const* y){
        detail::add192_192(x,y); // 384 + 192 = 192
    };
    
    template<>
    void add<unsigned long int,6>(unsigned long int * x,unsigned long int const y){
        detail::add384_64(x,&y); // 384 + 64 = 384
    }; 
    
    template<>
    void add<unsigned long int,3>(unsigned long int * x,unsigned long int const y){
        detail::add192_64(x,&y); // 384 + 64 = 192
    };

    //specialization sub
    template<>
    void sub<unsigned long int,6>(unsigned long int * x,unsigned long int const* y){
        detail::sub384_384(x,y); // 384 - 384 = 384
    }; 
    
    template<>
    void sub<unsigned long int,3>(unsigned long int * x,unsigned long int const* y){
        detail::sub192_192(x,y); // 384 - 192 = 192
    };
    
    template<>
    void sub<unsigned long int,6>(unsigned long int * x,unsigned long int const y){
        detail::sub384_64(x,&y); // 384 - 64 = 384
    }; 
    
    template<>
    void sub<unsigned long int,3>(unsigned long int * x,unsigned long int const y){
        detail::sub192_64(x,&y); // 384 - 64 = 192
    };
    
    //specialization mul    
    template<>
    void mul<unsigned long int,3>(unsigned long int * x,unsigned long int const y){
        detail::mul192_64(x,&y); // 384 * 64 = 384
    }; 
    
    template<>
    void mul<unsigned long int,6>(unsigned long int * x,unsigned long int const y){
        detail::mul384_64(x,&y); // 384 * 64 = 384
    }; 
        
    
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
