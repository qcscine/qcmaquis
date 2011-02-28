/*
 *  kernels.h
 *  Untitled
 *
 *  Created by Tim Ewart on 24.02.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#include "vli_matrix.h"
#include "definition.h"

namespace vli
{

	
/**
addition classic version on array 
*/
template <typename T>
void addition_classic_cpu(vli::vli_matrix<T>  const & x, vli::vli_matrix<T>  const &  y , vli::vli_matrix<T>  &  z)
{
	T carry_bit = 0;
	size_type num_rows = x.num_rows(); 
	size_type num_cols = x.num_columns();
	
	for (size_type i = 0 ; i < num_rows; i++) 
	{
		for(size_type j = 0 ; j < num_cols ; j++)
		{
			addition_kernel_cpu(&x(i,j), &y(i,j), &z(i,j));
			/* old version
			z(i,j)    += (x(i,j) + y(i,j));
			carry_bit  = z(i,j) >> LOG_BASE;
			z(i,j)    %= BASE;
			z(i+1,j)  += carry_bit; //MAYBE PB TO CHECK
			carry_bit = 0; //set to 0 for the nex columns
			*/
		}
	}
}
	

/**
addition classic version on array 
*/	
template <typename T>
void addition_kernel_cpu(T  const * x, T  const *  y , T * z)
{
	T carry_bit = 0;
	*z += *x + *y;
	carry_bit  = *z >> LOG_BASE;
	*z    %= BASE;
	*(z+1)+= carry_bit; //MAYBE PB TO CHECK
	carry_bit = 0; //set to 0 for the nex columns
}	
	
	
/**
addition second version
A. Avizienis. Signed-digit number representations for fast parallel arithmetic. IRE Transactions on electronic computers, vol 10, pages 389-400, 1961
No carry bit ^^'
*/	
template <typename T>
void addition_Avizienis_kernel_cpu(vli::vli_matrix<T>  const & x, vli::vli_matrix<T>  const &  y , vli::vli_matrix<T>  &  z)
{

	size_type num_rows = x.num_rows(); 
	size_type num_cols = x.num_columns();	
	
	
	T* t = (T*)malloc((num_rows+1)*sizeof(T)); //+1 for the step
	memset((void*)t,0,(num_rows+1)*sizeof(T));
	T* w = (T*)malloc(num_rows*sizeof(T));
	memset((void*)w,0,num_rows*sizeof(T));	
	T* s = (T*)malloc(num_rows*sizeof(T));
	memset((void*)s,0,num_rows*sizeof(T));

	for (size_type i = 0 ; i < num_rows; i++) 
	{
		#pragma omp parallel for firstprivate (i)
		for(size_type j = 0 ; j < num_cols ; j++)
		{
			s[i] = (x(i,j) + y(i,j));
/**
this three if could be optimize, to do ......
*/
			if(s[i] > BASE_MINUS2)
				t[i+1] = 1;
			if(s[i] < MINUS_BASE_PLUS2)
			    t[i+1] = -1;
			if(s[i]<=BASE_MINUS2 && s[i]>= MINUS_BASE_PLUS2)
				t[i+1] = 0;
		
			w[i] = s[i] - BASE*t[i+1];
		    z(i,j) = t[i] + w[i];
		}
	}
}
	
/**
multiplication classic version, efficiency O(n**2)
*/
template <typename T>
void multiplication_classic_kernel_cpu(vli::vli_matrix<T>  const & x, vli::vli_matrix<T>  const &  y , vli::vli_matrix<T>  &  z)	
{
	size_type num_rows = x.num_rows(); 
	size_type num_cols = x.num_columns();
	
	T* t1 = (T*)malloc((2*num_rows)*sizeof(T));
	memset((void*)t1,0,(2*num_rows)*sizeof(T));
	T* t2 = (T*)malloc((2*num_rows)*sizeof(T));
	memset((void*)t2,0,(2*num_rows)*sizeof(T));
	T* t3 = (T*)malloc((2*num_rows)*sizeof(T));
	memset((void*)t3,0,(2*num_rows)*sizeof(T));	
		
	T carry_bit = 0 ;
	size_type k,o;
	for (size_type i = 0 ; i < num_rows; i++) //loop on rows
	{
		k = i;
		o = i;
		for(size_type j = 0 ; j < num_cols ; j++) // loop on columns
		{
			
			for(; k < num_rows;k++) // loop on the local array for partial sum
			{
				t1[o]     = ((x(i,j) & MASK_DOWN)) * ((y(k,j)  & MASK_DOWN));
				t1[o+1]   = ((x(i,j) & MASK_DOWN)) * ((y(k,j)  & MASK_UP) >> LOG_BASE_HALF);
			
				t2[o+1]	  = ((x(i,j) & MASK_UP>> LOG_BASE_HALF ))   * ((y(k,j)  & MASK_DOWN));	
				t2[o+2]	  = ((x(i,j) & MASK_UP>> LOG_BASE_HALF ))   * ((y(k,j)  & MASK_UP) >> LOG_BASE_HALF);
			
				t3[o]	 += t2[o]+t1[o];
				carry_bit = t3[o] >> LOG_BASE;
				t3[o]    %= BASE;
				t3[k+1]    += carry_bit; 
				carry_bit = 0;
				
				o=o+2;
			}
		}
	}
}

	
/**
multiplication karatsuba version, efficiency O(n**1.58)
*/
template <typename T>
void multiplication_karatsuba_kernel_cpu(vli::vli_matrix<T>  const & x, vli::vli_matrix<T>  const &  y , vli::vli_matrix<T>  &  z)	
{
	// TO DO
}	

}