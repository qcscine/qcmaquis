/*
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 07.07.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */


#ifndef VLI_VECTOR_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP
#define VLI_VECTOR_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP

#include "detail/kernels_gpu.h"
#include "gpu/GpuManager.h"

namespace vli{

    template <class Polynomial> 
    class vector_polynomial_gpu;
    
    template <class BaseInt, int Size, int Order> 
    void inner_product_multiplication_gpu(
                           vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& a,
                           vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& b,
                           polynomial_gpu<vli_gpu<BaseInt, Size>, Order>& res)
    {
//        gpu::gpu_manager* GPU; //GPU is a singleton initialized into the main
        std::size_t NumThreads = 256 ; //GPU->instance().GetmaxThreadsPerBlock();
        
        vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> > inter;
        inter.resize(a.size()); //default size 1 so resize        
        detail::inner_product_vector(a.p(), b.p(),inter.p(),a.size(),NumThreads);
        detail::vector_reduction(inter.p(),res.p(),a.size());
    }
    
    
    
}// end namespace

#endif
