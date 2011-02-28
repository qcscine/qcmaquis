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
addition classic version 
*/
template <typename T>
void addition_classic_kernel_cpu(vli::vli_matrix<T>  const & x, vli::vli_matrix<T>  const &  y , vli::vli_matrix<T>  &  z)
{
	T carry_bit = 0;
	size_type num_rows = x.num_rows(); 
	size_type num_cols = x.num_columns();
	
	for (size_type i = 0 ; i < num_rows; i++) 
	{
		for(size_type j = 0 ; j < num_cols ; j++)
		{
			z(i,j)    += (x(i,j) + y(i,j));
			carry_bit  = z(i,j) >> LOG_BASE;
			z(i,j)    %= BASE;
			z(i+1,j)  += carry_bit; //MAYBE PB TO CHECK
			carry_bit = 0; //set to 0 for the nex columns
		}
	}
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

}