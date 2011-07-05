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

#include "common_macros.h"
#include <cassert>
#include <boost/static_assert.hpp>
#include <cstring>

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
	void addition_kernel_cpu(T* x, T const*  y)
	{
		T carry_bit = 0;
		*x += *y; 
		carry_bit  = *x >> LOG_BASE;
		*x    %= BASE; 
		*(x+1)+= carry_bit; //MAYBE PB TO CHECK
	}	
	
	/**
	 addition for the multiplication, just for a plus
	 */
	template <typename T>
	void addition_kernel2_cpu(T const* x, T const*  y , T* z)
	{
		T carry_bit = 0;
		*z = *x + *y; // <- no plus ><
		carry_bit  = *z >> LOG_BASE;
		*z    %= BASE;
		*(z+1)+= carry_bit; //MAYBE PB TO CHECK
	}
    
	/**
	 addition classic version on array 
	 */
	template <class T, int Size>
	void addition_classic_cpu(T* x,  T const*  y)
	{
        std::size_t size = Size;
		for(size_type i = 0; i < size ; ++i)
			addition_kernel_cpu((x+i), (y+i));
	}
	
	template <typename T>
	void multiplication_kernel_up_cpu(T const* x, T const*  y, T * r)	
	{
		*r	   = ((*x & MASK_UP) >> LOG_BASE_HALF ) * (*y & MASK_DOWN);	
		*(r+1) = ((*x & MASK_UP) >> LOG_BASE_HALF ) * ((*y & MASK_UP) >> LOG_BASE_HALF);
	}		
	
	template <typename T>
	void multiplication_kernel_down_cpu(T const* x, T const*  y, T * r)	
	{	
		*r     = (*x & MASK_DOWN) * (*y & MASK_DOWN);
		*(r+1) = (*x & MASK_DOWN) * ((*y & MASK_UP) >> LOG_BASE_HALF);
	}
	
	template <typename T>
	void multiplication_kernel_base_reshaping_cpu(T const * a, T  const *  b, T * r)	
	{			
		int q1,q2;
		int r1,r2;
		q1 = q2 = r1 =r2 = 0;
		
		q1 = (*(a+1) + *b)/BASE_HALF;
		r1 = (*(a+1) + *b)%BASE_HALF;
		r1 = r1 * BASE_HALF;
		q2 = (r1 + *a)/BASE; 
		r2 = (r1 + *a)%BASE;
		*r  = r2;
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
		 = (q1+q2 + Xl*Yl)*BASE + r2  (multiplication_kernel_base_reshaping)
		 */
		multiplication_kernel_down_cpu(x,y, a);
		multiplication_kernel_up_cpu(x,y, b);
		multiplication_kernel_base_reshaping_cpu(a,b, r);
	}		
	
	/**
	 multiplication classic version, efficiency O(n**2)
	 */
	template <typename BaseInt, int Size>
	void multiplication_classic_cpu(BaseInt* x, const BaseInt* y)	
	{
        std::size_t size = Size;

		BaseInt r[2] = {0,0};	//for local block calculation
		BaseInt m = 0;
		
        vli::vli_cpu<BaseInt, Size> inter; 
		
		for (size_type i = 0 ; i < size; i++)
		{
			for(size_type k = 0 ; k < size ; k++) // loop on numbers for multiplication the classical multiplication
			{	
				m = static_cast<BaseInt>(k + i);
				multiplication_block_cpu( &x[i], &y[k], &(r[0]));
				addition_kernel_cpu(&inter[m],&r[0]);
				addition_kernel_cpu(&inter[m+1],&r[1]);
			}
		}
		
        memcpy((void*)x,(void*)&inter[0],Size*sizeof(BaseInt));
	}
    
	/***************************** choses exotiques *************************************/
	
	/**
	 multiplication karatsuba version, efficiency O(n**1.58)
	 */
    //	template <typename T>
    //	void multiplication_karatsuba_kernel_cpu(vli::vli_matrix<T>  const & x, vli::vli_matrix<T>  const &  y , vli::vli_matrix<T>  &  z)	
    //	{
    //		// TO DO
    //	}
	
    /**
	 addition second version
	 A. Avizienis. Signed-digit number representations for fast parallel arithmetic. IRE Transactions on electronic computers, vol 10, pages 389-400, 1961
	 No carry bit ^^'
	 */	
    //	template <typename T>
    //	void addition_Avizienis_kernel_cpu(vli::vli_matrix<T>  const & x, vli::vli_matrix<T>  const &  y , vli::vli_matrix<T>  &  z)
    //	{
    //		
    //		size_type num_rows = x.num_rows(); 
    //		size_type num_cols = x.num_columns();	
    //		
    //		/**
    //		 TO DO : REMOVE ARRAY AND ISOLATE THE KERNEL 
    //		 */
    //		T* t = (T*)malloc((num_rows+1)*sizeof(T)); //+1 for the step
    //		memset((void*)t,0,(num_rows+1)*sizeof(T));
    //		T* w = (T*)malloc(num_rows*sizeof(T));
    //		memset((void*)w,0,num_rows*sizeof(T));	
    //		T* s = (T*)malloc(num_rows*sizeof(T));
    //		memset((void*)s,0,num_rows*sizeof(T));
    //		
    //		for (size_type i = 0 ; i < num_rows; i++) 
    //		{
    //#pragma omp parallel for firstprivate (i)
    //			for(size_type j = 0 ; j < num_cols ; j++)
    //			{
    //				s[i] = (x(i,j) + y(i,j));
    //				/**
    //				 this three if could be optimize, to do ...... il else statement
    //				 */
    //				if(s[i] > BASE_MINUS2)
    //					t[i+1] = 1;
    //				if(s[i] < MINUS_BASE_PLUS2)
    //					t[i+1] = -1;
    //				if(s[i]<=BASE_MINUS2 && s[i]>= MINUS_BASE_PLUS2)
    //					t[i+1] = 0;
    //				
    //				w[i] = s[i] - BASE*t[i+1];
    //				z(i,j) = t[i] + w[i];
    //			}
    //		}
    //	}
}

#endif //VLI_KERNELS_CPU_HPP
