/*
 *  kernel_SIMD.h
 *  dmrg
 *
 *  Created by Tim Ewart on 22/01/11.
 *  Copyright 2011 Université de Genève. All rights reserved.
 *
 */

#ifndef KERNEL_SIMD_H
#define KERNEL_SIMD_H

#include <emmintrin.h>

/**
	Only for double run on INTEL/AMD with gcc 4.4.4 write with SSE2 to be compatible
	with a large set of precessors
*/

void kernel_vector_scalar(double * p1, double const * p2, double  scalaire, std::size_t n )
{
	/**
	 vectors register 
	*/
	 
	__m128d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8;
	std::size_t m;
	
	double scalaire_array[2] __attribute__ ((aligned (16))) = {scalaire, scalaire};
	
	xmm0 = _mm_load_pd(&scalaire_array[0]);
	
	if (n/8 >= 1)
	{
		
		m  =  8*std::size_t(n/8);
		
		for(std::size_t i = 0; i < m ; i = i+8)
		{
			std::cout << "simd" << std::endl;
			/**
			 objective :  *(p1++) += *(p2++) * wblock_t;
			 unroll + SSE2
			 */
			//load p1
			xmm1 = _mm_load_pd(p1+i);
			xmm2 = _mm_load_pd(p1+(i+2));
			xmm3 = _mm_load_pd(p1+(i+4));
			xmm4 = _mm_load_pd(p1+(i+6));
			//load p2
			xmm5 = _mm_load_pd(p2+i);
			xmm6 = _mm_load_pd(p2+(i+2));
			xmm7 = _mm_load_pd(p2+(i+4));
			xmm8 = _mm_load_pd(p2+(i+6));
			//*(p1++) += *(p2++) * wblock_t;
			xmm5 = _mm_mul_pd(xmm5, xmm0);
			xmm6 = _mm_mul_pd(xmm6, xmm0);
			xmm7 = _mm_mul_pd(xmm7, xmm0);
			xmm8 = _mm_mul_pd(xmm8, xmm0);
			
			xmm1 = _mm_add_pd(xmm1, xmm5);
			xmm2 = _mm_add_pd(xmm2, xmm6);
			xmm3 = _mm_add_pd(xmm3, xmm7);
			xmm4 = _mm_add_pd(xmm4, xmm8);
			//set p1 et p2 
			_mm_store_pd(p1+i  , xmm1);
			_mm_store_pd(p1+(i+2), xmm2);			
			_mm_store_pd(p1+(i+4), xmm3);
			_mm_store_pd(p1+(i+6), xmm4);
			
		}
		
	}
	
	//finish classic
	for (std::size_t i = m; i < n; i++)
	{
		p1[i] += p2[i]*scalaire;
	}
	
};




#endif