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

    template<typename BaseInt, std::size_t Size1, std::size_t Size2>
    void add(BaseInt * x, BaseInt const* y, BaseInt const* z); 

    //substraction
    template<typename BaseInt, std::size_t Size>
    void sub(BaseInt * x, BaseInt const b);
    
    template<typename BaseInt, std::size_t Size>
    void sub(BaseInt * x, BaseInt const* y); 

    //multiplication
    template<typename BaseInt, std::size_t Size>
    void mul(BaseInt * x,BaseInt const b);

    template<typename BaseInt, std::size_t Size>
    void mul(BaseInt * x,BaseInt const* y);

    //????_assign functions
    template <class BaseInt, std::size_t Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, vli_cpu<BaseInt,Size> const& vli_b ){
        add<BaseInt,Size>(&vli_a[0],&vli_b[0]);
    }
     
    template <class BaseInt, std::size_t Size>
    void plus_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b ){  
        add<BaseInt,Size>(&vli_a[0],b);
    }
    
    template <class BaseInt, std::size_t Size1, std::size_t Size2>
    void addition_extension(vli_cpu<BaseInt,Size2> & vli_a, vli_cpu<BaseInt,Size1> const& vli_b, vli_cpu<BaseInt,Size1> const& vli_c){
        add<BaseInt, Size1, Size2>(&vli_a[0],&vli_b[0],&vli_c[0]);
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
    void multiplies_assign( vli_cpu<BaseInt, Size>& vli_a , vli_cpu<BaseInt,Size> const & vli_b){ 
        mul<BaseInt,Size>(&vli_a[0],&vli_b[0]);
    }

    template <class BaseInt, std::size_t Size>
    void multiplies_assign(vli_cpu<BaseInt,Size> & vli_a, BaseInt const b){
        mul<BaseInt,Size>(&vli_a[0],b);
    }

    template <class BaseInt, std::size_t Size>
    void multiplies(vli_cpu<BaseInt, 2*Size>& vli_res , vli_cpu<BaseInt,Size> const & vli_a, vli_cpu<BaseInt,Size> const & vli_b){
        detail::mul384_192_192(&vli_res[0],&vli_a[0],&vli_b[0]);
    }
    
    
    template <class BaseInt, std::size_t Size>
    void multiply_add_assign(vli_cpu<BaseInt, 2*Size>& vli_res , vli_cpu<BaseInt,Size> const & vli_a, vli_cpu<BaseInt,Size> const & vli_b){
        detail::muladd384_192_192(&vli_res[0],&vli_a[0],&vli_b[0]);
    }

    /* ---------------------------------------------------- Begin Addition specialization ---------------------------------------------------- */

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
        void add<unsigned long int,BOOST_PP_ADD(n,2)>(unsigned long int* x,unsigned long int const y){ \
        detail::NAME_ADD_NBITS_PLUS_64BITS(n)(x,&y); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_add_nbits_64bits, ~)
    #undef FUNCTION_add_nbits_64bits
    // specialization extention addition 
    #define FUNCTION_add_nbits_nminus1bits(z, n, unused) \
        template<> \
        void add<unsigned long int,BOOST_PP_ADD(n,1),BOOST_PP_ADD(n,2)>(unsigned long int* x,unsigned long int const* y, unsigned long int const* w){ \
        detail::NAME_ADD_NBITS_PLUS_NMINUS1BITS(n)(x,y,w); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION_MINUS_ONE, FUNCTION_add_nbits_nminus1bits, ~)
    #undef FUNCTION_add_nbits_mninus1bits
    /* ---------------------------------------------------- End Addition specialization ---------------------------------------------------- */

    /* ---------------------------------------------------- Begin Substraction specialization ---------------------------------------------------- */
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
        void sub<unsigned long int,BOOST_PP_ADD(n,2)>(unsigned long int* x,unsigned long int const y){ \
        detail::NAME_SUB_NBITS_MINUS_64BITS(n)(x,&y); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_sub_nbits_64bits, ~)
    #undef FUNCTION_sub_nbits_64bits

    /* ---------------------------------------------------- end Substraction specialization ---------------------------------------------------- */

    /* ---------------------------------------------------- Begin Multipliation specialization ---------------------------------------------------- */
    //specialization mul    
    //specialization mulnbits_64bits, until 512 bits
    #define FUNCTION_mul_nbits_64bits(z, n, unused) \
        template<> \
        void mul<unsigned long int,BOOST_PP_ADD(n,2)>(unsigned long int* x,unsigned long int const y){ \
        detail::NAME_MUL_NBITS_64BITS(n)(x,&y); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_mul_nbits_64bits, ~)
    #undef FUNCTION_mul_nbits_64bits
    //specialization mulnbits_nbits, until 512 bits
    #define FUNCTION_mul_nbits_nbits(z, n, unused) \
        template<> \
        void mul<unsigned long int,BOOST_PP_ADD(n,2)>(unsigned long int* x,unsigned long int const* y){ \
        detail::NAME_MUL_NBITS_NBITS(n)(x,y); \
        }; \

    BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_mul_nbits_nbits, ~)
    #undef FUNCTION_mul_nbits_nbits
    /* ---------------------------------------------------- end Multiplicatio specialization ---------------------------------------------------- */

    
} //namespace vli

#endif //VLI_NUMBER_CPU_FUNCTION_HOOKS_HPP
