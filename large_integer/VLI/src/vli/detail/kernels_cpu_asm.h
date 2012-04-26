//
//  kernels_cpu_asm.h
//  VLI_ASM
//
//  Created by Tim Ewart on 22.03.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef VLI_ASM_KERNELS_CPU_ASM_H
#define VLI_ASM_KERNELS_CPU_ASM_H 

#include "vli/utils/macro.h" 

namespace vli{
    namespace detail{

        // C first number output #bits, second and third input #bits
        //Addition
        // new functions type : VLI<n*64> + VLI<n*64> : add128_128, add192_192 ...
        #define FUNCTION_add_nbits_nbits(z, n, unused) \
            void NAME_ADD_NBITS_PLUS_NBITS(n)(unsigned long int* x, unsigned long int const* y);

        BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_add_nbits_nbits, ~)
        #undef FUNCTION_add_nbits_nbits
        /* ------------------------------------------------------- */
        //new functions type : VLI<n*64> + VLI<64> : add192_64 ...
        #define FUNCTION_add_nbits_64bits(z, n, unused) \
            void NAME_ADD_NBITS_PLUS_64BITS(n)(unsigned long int* x, unsigned long int const* y);

        BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_add_nbits_64bits, ~)
        #undef FUNCTION_add_nbits_64bits
        /* ------------------------------------------------------- */
        //new functions type : VLI<n*64> + VLI<64> : add128_64, add128_64 ...
        #define FUNCTION_add_nbits_nminus1bits(z, n, unused) \
            void NAME_ADD_NBITS_PLUS_NMINUS1BITS(n)(unsigned long int* x, unsigned long int const* y);

        BOOST_PP_REPEAT(MAX_ITERATION_MINUS_ONE, FUNCTION_add_nbits_nminus1bits, ~)
        #undef FUNCTION_add_nbits_nminus1bits
        /* ------------------------------------------------------- */
        //substraction
        #define FUNCTION_sub_nbits_nbits(z, n, unused) \
            void NAME_SUB_NBITS_MINUS_NBITS(n)(unsigned long int* x, unsigned long int const* y);

        BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_sub_nbits_nbits, ~)
        #undef FUNCTION_sub_nbits_nbits
        /* ------------------------------------------------------- */
        #define FUNCTION_sub_nbits_64bits(z, n, unused) \
            void NAME_SUB_NBITS_MINUS_64BITS(n)(unsigned long int* x, unsigned long int const* y);

        BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_sub_nbits_64bits, ~)
        #undef FUNCTION_sub_nbits_64bits
        /* ------------------------------------------------------- */
        #define FUNCTION_sub_nbits_nminus1bits(z, n, unused) \
            void NAME_SUB_NBITS_MINUS_NMINUS1BITS(n)(unsigned long int* x, unsigned long int const* y);

        BOOST_PP_REPEAT(MAX_ITERATION_MINUS_ONE, FUNCTION_sub_nbits_nminus1bits, ~)
        #undef FUNCTION_sub_nbits_nminus1bits
        /* ------------------------------------------------------- */
        //multiplication
        #define FUNCTION_mul_nbits_64bits(z, n, unused) \
            void NAME_MUL_NBITS_64BITS(n)(unsigned long int* x,unsigned long int const* y);

        BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_mul_nbits_64bits, ~)
        #undef FUNCTION_mul_nbits_64bits
        /* ------------------------------------------------------- */
        #define FUNCTION_mul_nbits_nbits(z, n, unused) \
            void NAME_MUL_NBITS_NBITS(n)(unsigned long int* x,unsigned long int const* y);

        BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_mul_nbits_nbits, ~)
        #undef FUNCTION_mul_nbits_nbits

        void mul384_192_192(unsigned long int * x,unsigned long int const* y,unsigned long int const* z); // 192*192 = 384
        void muladd384_192_192(unsigned long int * x,unsigned long int const* y,unsigned long int const* z); // 384 += 192*192
    }
}
        
#endif
