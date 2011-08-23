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
#include <iostream>
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
	void multiplication_block_cpu(T const* x, T  const*  y, T * r )	
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
		kernels::kernels_multiplication_block_down(x,y, a);
        kernels::kernels_multiplication_block_up(x,y, b);
        kernels::kernels_multiplication_base_reshaping(a,b,r);
	}		
	
	/**
    Only used for the polynomial inner product
	 multiplication classic version, efficiency O(n**2)
     the size of res is 2*size
	 */
	template <typename BaseInt, std::size_t Size>
	void multiplication_classic_cpu(BaseInt * res, BaseInt const* x, BaseInt const* y)	
	{
		BaseInt r[2] = {0,0};	//for local block calculation
        std::size_t m(0);
        
		for(std::size_t i = 0 ; i < Size; ++i)
		{
			for(std::size_t k = 0 ; k < Size; ++k) // loop on numbers for multiplication the classical multiplication
			{	
                m = k + i;
				multiplication_block_cpu( &x[i], &y[k], &(r[0]));
				addition_kernel_cpu(&res[m],&r[0]);
				addition_kernel_cpu(&res[m+1],&r[1]);
			}
		}
    }
    
    
    /** 
    truncated multiplication
    */
    template <typename BaseInt, std::size_t Size>
	void multiplication_classic_truncate_cpu(BaseInt * res, BaseInt const* x, BaseInt const* y)	
	{
       kernels::kernels_multiplication_classic_truncate<BaseInt,Size>(res,x,y);
    }
    
    /** 
     This multiplication is done from right to left
    */
    template <typename BaseInt, std::size_t Size>
	void multiplication_classic_cpu_number(BaseInt* x, BaseInt a)	
	{      
        kernels::kernels_multiplication_number<BaseInt,Size>(x,a);  
    }

}

#endif //VLI_KERNELS_CPU_HPP
