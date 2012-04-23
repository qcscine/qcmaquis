#ifndef VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
#define VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP

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

    //????_assign functions
    template <class BaseInt, std::size_t Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b ){
        add<BaseInt,Size>(&vli_a[0],&vli_b[0]);
    }
     
    template <class BaseInt, std::size_t Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b ){  
        add<BaseInt,Size>(&vli_a[0],b);
    }
    
    template <class BaseInt, std::size_t Size>
    void minus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b ){
        sub<BaseInt,Size>(&vli_a[0],&vli_b[0]);
    }
    
    template <class BaseInt, std::size_t Size>
    void minus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b ){  
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
    
    template <class BaseInt, std::size_t Size>
    void multiply_add_assign(vli_cpu<BaseInt, 2*Size>& vli_res , vli_cpu<BaseInt,Size> const & vli_a, vli_cpu<BaseInt,Size> const & vli_b){
        detail::muladd384_192_192(&vli_res[0],&vli_a[0],&vli_b[0]);
    }

    //specialization addnbits_nbits, until 512 bits
    #define FUNCTION_add_nbits_nbits(z, n, unused) \
        template<> \
        void add<unsigned long int,BOOST_PP_ADD(n,2)>(unsigned long int* x,unsigned long int const* y){ \
        detail::NAME_ADD_NBITS_PLUS_NBITS(n)(x,y); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_add_nbits_nbits, ~)
    #undef FUNCTION_add_nbits_nbits

    //specialization addnbits_64bits, until 512 bits
    #define FUNCTION_add_nbits_64bits(z, n, unused) \
        template<> \
        void add<unsigned long int,BOOST_PP_ADD(n,3)>(unsigned long int* x,unsigned long int const y){ \
        detail::NAME_ADD_NBITS_PLUS_64BITS(n)(x,&y); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION_MINUS_ONE, FUNCTION_add_nbits_64bits, ~)
    #undef FUNCTION_add_nbits_64bits

    //specialization subnbits_nbits, until 512 bits

    #define FUNCTION_sub_nbits_nbits(z, n, unused) \
        template<> \
        void sub<unsigned long int,BOOST_PP_ADD(n,2)>(unsigned long int* x,unsigned long int const* y){ \
        detail::NAME_SUB_NBITS_MINUS_NBITS(n)(x,y); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_sub_nbits_nbits, ~)
    #undef FUNCTION_sub_nbits_nbits

    //specialization subnbits_64bits, until 512 bits
    #define FUNCTION_sub_nbits_64bits(z, n, unused) \
        template<> \
        void sub<unsigned long int,BOOST_PP_ADD(n,3)>(unsigned long int* x,unsigned long int const y){ \
        detail::NAME_SUB_NBITS_MINUS_64BITS(n)(x,&y); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION_MINUS_ONE, FUNCTION_sub_nbits_64bits, ~)
    #undef FUNCTION_sub_nbits_64bits
    //specialization mul    
    template<>
    void mul<unsigned long int,3>(unsigned long int * x,unsigned long int const y){
        detail::mul192_64(x,&y); // 384 * 64 = 384
    }; 
    
    void totoa(unsigned long int * x,unsigned long int const y){
        detail::mul192_64(x,&y); // 384 * 64 = 384        
    }    
    void totob(unsigned long int * x,unsigned long int const y){
        detail::mul192_64b(x,&y); // 384 * 64 = 384        
    }
    
    template<>
    void mul<unsigned long int,6>(unsigned long int * x,unsigned long int const y){
        detail::mul384_64(x,&y); // 384 * 64 = 384
    }; 
        
    
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
