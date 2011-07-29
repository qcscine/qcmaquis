/*
 *  monome_gpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 07.07.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */


#ifndef VLI_VECTOR_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP
#define VLI_VECTOR_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP

#include "gpu/GpuManager.h"
#include "detail/kernels_gpu.h"

namespace vli{
    

    template <class Polynomial> 
    class vector_polynomial_gpu;
    
    template <class BaseInt, int Size, int Order> 
    void inner_product_multiplication_gpu(
                           vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& a,
                           vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >  const& b,
                           vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> >   & res)
    {
        //As inter is dynamic I prefer make the alloc by the CPU ....
        vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> > inter;
        inter.resize(a.size());
        detail::inner_product_vector_gpu(a.p(), b.p(),inter.p(), Size, Order, a.size());
        detail::vector_reduction_gpu(inter.p(),res.p(),Size, Order, inter.size());
    }
    
    
    
}// end namespace

#endif
