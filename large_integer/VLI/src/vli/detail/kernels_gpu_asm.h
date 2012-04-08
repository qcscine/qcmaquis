//
//  kernels_gpu_asm.h
//  vli
//
//  Created by Timothée Ewart on 08/04/12.
//  Copyright (c) 2012 IBM. All rights reserved.
//

//#ifndef vli_kernels_gpu_asm_h
//#define vli_kernels_gpu_asm_h

namespace vli{
    namespace detail{
        __device__ void add192_192(unsigned long int* x, unsigned long int const* y);
    }
}

//#endif
