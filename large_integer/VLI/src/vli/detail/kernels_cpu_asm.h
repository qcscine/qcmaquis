//
//  kernels_cpu_asm.h
//  VLI_ASM
//
//  Created by Tim Ewart on 22.03.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

//#ifndef VLI_ASM_KERNELS_CPU_ASM_H
//#define VLI_ASM_KERNELS_CPU_ASM_H 

#include <iostream>

namespace vli{
    namespace detail{

        // C first number output #bits, second and third input #bits
        //Addition
        void add384_64(unsigned long int * x,unsigned long int const* y);            
        void add384_384(unsigned long int * x,unsigned long int const* y);            
        void add192_192(unsigned long int * x,unsigned long int const* y);            
        void add192_64(unsigned long int * x,unsigned long int const* y);            
        
        //substraction
        void sub384_64(unsigned long int * x,unsigned long int const* y);            
        void sub384_384(unsigned long int * x,unsigned long int const* y);            
        void sub192_192(unsigned long int * x,unsigned long int const* y);            
        void sub192_64(unsigned long int * x,unsigned long int const* y);            

        //multiplication
        void mul384_64(unsigned long int * x,unsigned long int const* y); // 384 * 64 = 384
        void mul384_192_192(unsigned long int * x,unsigned long int const* y,unsigned long int const* z); // 192*192 = 384
        void mul192_192(unsigned long int * x,unsigned long int const* y); // 192 * 192 = 192
        void mul192_64(unsigned long int * x,unsigned long int const* y); // 192 * 64 = 192
    }
}
        
//#endif
