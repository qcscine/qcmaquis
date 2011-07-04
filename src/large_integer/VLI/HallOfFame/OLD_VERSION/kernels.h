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
}	

/**
addition for the multiplication, just for a plus
*/
template <typename T>
void addition_kernel2_cpu(T  const * x, T  const *  y , T * z)
{
	T carry_bit = 0;
	*z = *x + *y; // <- no plus ><
	carry_bit  = *z >> LOG_BASE;
	*z    %= BASE;
	*(z+1)+= carry_bit; //MAYBE PB TO CHECK
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
	
	/**
	 TO DO : REMOVE ARRAY AND ISOLATE THE KERNEL 
	 */
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
this three if could be optimize, to do ...... il else statement
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

	
template <typename T>
void multiplication_kernel_up_cpu(T  const * x, T  const *  y, T * r)	
{
	*r	    = ((*x & MASK_UP) >> LOG_BASE_HALF ) * (*y & MASK_DOWN);	
	*(r+1)	= ((*x & MASK_UP) >> LOG_BASE_HALF ) * ((*y & MASK_UP) >> LOG_BASE_HALF);
}		
	
template <typename T>
void multiplication_kernel_down_cpu(T  const * x, T  const *  y, T * r)	
{	
	*r     = (*x & MASK_DOWN) * (*y & MASK_DOWN);
	*(r+1) = (*x & MASK_DOWN) * ((*y & MASK_UP) >> LOG_BASE_HALF);
}

template <typename T>
void multiplication_kernel_base_reshaping_cpu(T  const * a, T  const *  b, T * r)	
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
void multiplication_block_cpu(T  const * x, T  const *  y, T * r)	
{
	int a[2] = {0,0};
	int b[2] = {0,0};
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
template <typename T>
void multiplication_classic_kernel_cpu(vli::vli_matrix<T>  const & x, vli::vli_matrix<T>  const &  y , vli::vli_matrix<T>  &  z)	
{
	size_type num_rows = x.num_rows(); 
	size_type num_cols = x.num_columns();
	int r[2] = {0,0};	//for local block calculation
	int m =0;
	
	for (size_type i = 0 ; i < num_rows; i++) //loop on rows
	{
		for(size_type j = 0 ; j < num_cols ; j++) // loop on columns
		{	
			for(size_type k = 0 ; k < num_rows ; k++) // loop on number for multiplication the classical multiplication
			{	
				m = k + i;
				multiplication_block_cpu(&x(i,j), &y(k,j), r);
				addition_kernel2_cpu(&z(m,j),&r[0],&z(m,j));
				addition_kernel2_cpu(&z(m+1,j),&r[1],&z(m+1,j));
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