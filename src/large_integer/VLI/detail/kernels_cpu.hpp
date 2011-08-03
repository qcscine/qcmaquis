/*
 *  kernels_cpu.hpp
 *  Untitled
 *
 *  Created by Tim Ewart on 24.02.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#ifndef VLI_KERNELS_CPU_HPP 
#define VLI_KERNELS_CPU_HPP

#include <cassert>
#include <boost/static_assert.hpp>
#include <cstring>
#include "detail/bit_masks.hpp"

namespace vli
{
    /*
     * Wellcome in the pointer jungle
     */
    template <typename BaseInt, int Size>
    class vli_cpu;    
	
    /**
	 addition classic version on array 
	 */	
	template <typename T>
	inline void addition_kernel_cpu(T* x, T const*  y)
	{
		*x += *y; 
		*(x+1) += *x >> data_bits<T>::value; //carry bit
		*x    &= data_mask<T>::value; // Remove the carry bit
	}	
    
	/**
	 addition classic version on array 
	 */
	template <class T, std::size_t Size>
	void addition_classic_cpu(T* x,  T const*  y)
	{
		for(std::size_t i = 0; i < Size-1 ; ++i)
			addition_kernel_cpu((x+i), (y+i));
		*(x+Size-1) += *(y+Size-1);
        *(x+Size-1) &= base<T>::value + data_mask<T>::value;
	}
    
	template <typename T>
	void multiplication_kernel_up_cpu(T const* x, T const*  y, T * r)	
	{
		*r	   = ((*x & mask_up<T>::value) >> (data_bits<T>::value/2) ) * (*y & mask_down<T>::value);	
		*(r+1) = ((*x & mask_up<T>::value) >> (data_bits<T>::value/2) ) * ((*y & mask_up<T>::value) >> (data_bits<T>::value/2));
	}		
	
	template <typename T>
	void multiplication_kernel_down_cpu(T const* x, T const*  y, T * r)	
	{	
		*r     = (*x & mask_down<T>::value) * (*y & mask_down<T>::value);
		*(r+1) = (*x & mask_down<T>::value) * ((*y & mask_up<T>::value) >> (data_bits<T>::value/2));
	}
	
	template <typename T>
	void multiplication_kernel_base_reshaping_cpu(T const * a, T  const *  b, T * r)	
	{			
		T q1 = (*(a+1) + *b) >> (data_bits<T>::value/2);
		T r1 = (*(a+1) + *b) & mask_down<T>::value;
		r1 *= base_half<T>::value;
		T q2 = (r1 + *a) >> data_bits<T>::value; 
		T r2 = (r1 + *a) & data_mask<T>::value;
		*r = r2;
		*(r+1) = q1 + q2 + *(b+1);
	}
	
	template <typename T>
	void multiplication_block_cpu(T* x, T  const*  y, T * r )	
	{
		T a[2] = {0,0};
		T b[2] = {0,0};
		/**
		 Divide and conquer algo (see my notes - for the euclidian division tips)
		 X <=> Xl Xr (half of the binary number)
		 x  Y <=> Yl Yr (half of the binary number)
		 -------
		 = 2^n XlYl + 2^(n/2) (XlYr + XrYl) + XrYr (multiplication_kernel_cpu_down and multiplication_kernel_cpu_up)
		 ------- 
		 = (q1+q2 + Xl*Yl)*base<T>::value + r2  (multiplication_kernel_base_reshaping)
		 */
		multiplication_kernel_down_cpu(x,y, a);
		multiplication_kernel_up_cpu(x,y, b);
		multiplication_kernel_base_reshaping_cpu(a,b, r);
	}		
	
	/**
	 multiplication classic version, efficiency O(n**2)
	 */
	template <typename BaseInt, std::size_t Size>
	void multiplication_classic_cpu(BaseInt* x, const BaseInt* y)	
	{
		BaseInt r[2] = {0,0};	//for local block calculation
        BaseInt inter[2*Size+1]; // largeur for avoid overflow
        std::size_t DynamicSize = Size;
      
        //for unroll
        for(std::size_t i=0; i<(2*Size+1); ++i)
            inter[i]=0;
        
		for(std::size_t i = 0 ; i < Size; ++i)
		{
            /** 
             WE do not calculate all digits so DANGEROUS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
             */
			for(std::size_t k = 0 ; k < DynamicSize ; ++k) // loop on numbers for multiplication the classical multiplication
			{	
                std::size_t m = k + i;
				multiplication_block_cpu( &x[i], &y[k], &(r[0]));
				addition_kernel_cpu(&inter[m],&r[0]);
				addition_kernel_cpu(&inter[m+1],&r[1]);        
			}
            DynamicSize --;
		}
        for(std::size_t i=0; i<Size; ++i) // we keep only what we need ... think how to remove it
            x[i] = inter[i];
	}
}

#endif //VLI_KERNELS_CPU_HPP
