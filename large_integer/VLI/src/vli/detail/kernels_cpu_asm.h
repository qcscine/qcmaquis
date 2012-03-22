//
//  kernels_cpu_asm.h
//  VLI_ASM
//
//  Created by Tim Ewart on 22.03.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef VLI_ASM_KERNELS_CPU_ASM_H
#define VLI_ASM_KERNELS_CPU_ASM_H 


namespace vli{
    namespace detail{
//Addition
    void add192(unsigned long int * x,unsigned long int const* y);            
    void add64(unsigned long int * x,unsigned long int const* y);            

//substraction
    void sub192(unsigned long int * x,unsigned long int const* y);        
    void sub64(unsigned long int * x,unsigned long int const* y);        

//multiplication
    void mul192(unsigned long int * x,unsigned long int const* y,unsigned long int const* z);
    void mul192t(unsigned long int * x,unsigned long int const* y);
    void mul64(unsigned long int * x,unsigned long int const* y);    
    }
}
        
#endif
