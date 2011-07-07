/*
 *  monome_gpu.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 07.07.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */


#ifndef VLI_VECTOR_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP
#define VLI_VECTOR_POLYNOMIAL_GPU_FUNCTION_HOOKS_HPP

#include "detail/kernels_gpu.h"

namespace vli{
    

    template <class Polynomial> 
    class vector_polynomial_gpu;
    
    template <class BaseInt, int Size, int Order>
    vector_polynomial_gpu<polynomial_gpu<vli_gpu<BaseInt, Size>, Order> > inner_product_gpu(
                                                                                    vector_polynomial_gpu<polynomial<vli_gpu<BaseInt, Size>, Order> >  const& a,
                                                                                    vector_polynomial_gpu<polynomial<vli_gpu<BaseInt, Size>, Order> >  const& b,
                                                                                    vector_polynomial_gpu<polynomial<vli_gpu<BaseInt, Size>, Order> >   & res)
    {
      //  detail::inner_product_vector_gpu(a.p(), b.p(), res.p(), Size, Order, a.size()); 
    }
    
    
    
    
}// end namespace

#endif