//
//  kernels_cpu_gpu.hpp
//  vli
//
//  Created by Timothée Ewart on 22/08/11.
//  Copyright 2011 University of Geneva. All rights reserved.
//

#ifndef VLI_KERNELS_CPU_GPU_HPP 
#define VLI_KERNELS_CPU_GPU_HPP 

#include "detail/bit_masks.hpp"


/**
*
* All these functions are kernels, they are compiled by CPU and GPU
* take care of the compatibility !
*/


namespace vli
{
    /*
     * Wellcome in the pointer jungle
     */

    namespace detail
    {
    /**
    *
    * BEGINING ADDITION ALGO BETWEEN TWO VLI
    *
    * As in your childhood, we propagate a carry bit
    */
    
    /**
       addition classic version on array  (VLI)
    */
	template <class T, std::size_t Size>
	void kernels_addition_classic(T* x,  T const*  y)
	{
		for(std::size_t i = 0; i < Size-1 ; ++i)
			kernels_addition_block((x+i), (y+i));
     
		*(x+Size-1) += *(y+Size-1);
        *(x+Size-1) &= base<T>::value + data_mask<T>::value;
	}
    
    /**
	 addition classic on a block 
	 */	
    template <typename T>
	inline void kernels_addition_block(T* x, T const*  y)
	{
		*x     += *y; 
		*(x+1) += *x >> data_bits<T>::value; //carry bit
		*x     &= data_mask<T>::value; // Remove the carry bit
	}	
    /**
    * END ADDITION BETWEEN TWO VLI
    */

    /**
    *
    * BEGINING MULTIPLICATION ALGO BETWEEN TWO VLI
    * 
    * vli multiplication inner functions, main idea devide/conqueur algo
    * 2 steps multiplies half of the type (e.g. int  - long int)
    * by the kernels_multiplication_block_up/down functions
    * make the final addition between block taking care og the local base
    * kernels_multiplication_base_reshaping
    *
    * the nultiplication is truncated thus we save time
    * and int * int = int not int * int = long int !!
    */
    template <typename BaseInt, std::size_t Size>
    void kernels_multiplication_classic_truncate(BaseInt * res, BaseInt const* x, BaseInt const* y)	
    {
    /**
    * As your childhood, but DC algo when you mutiply block to avoir to lose data
    */
        BaseInt r[2] = {0,0};	//for local block calculation
        std::size_t m(0);
        std::size_t DynamicSize(Size);

        for(std::size_t i = 0 ; i < Size; ++i)
        {
            for(std::size_t k = 0 ; k < DynamicSize; ++k) // loop on numbers for multiplication the classical multiplication
            {
                m = k + i;
                kernels_multiplication_block( &x[i], &y[k], &(r[0]));
                kernels_addition_block(&res[m],&r[0]);
                kernels_addition_block(&res[m+1],&r[1]);
            }
            DynamicSize--;
         }
    }

    template <typename T>
    inline void kernels_multiplication_block(T const* x, T const * y, T* r)
    {
        T a[2] = {0,0};
        T b[2] = {0,0};
        /**
         Divide and conquer algo (see my notes - for the euclidian division tips)
          X  <=> Xl Xr (half of the binary number)
          Y  <=> Yl Yr (half of the binary number)
         -------
          = 2^n XlYl + 2^(n/2) (XlYr + XrYl) + XrYr (kernels_multiplication_block_down and kernels_multiplication_block_up)
         ------- 
          = (q1+q2 + Xl*Yl)*BASE + r2  (multiplication_kernel_base_reshaping)
        */
        kernels_multiplication_block_down(x,y,a);
        kernels_multiplication_block_up(x,y,b);
        kernels_multiplication_base_reshaping(a,b, r);
    }
 
    template <typename T>
    inline void kernels_multiplication_block_down(T const* x, T const*  y, T * r)	
    {
        *r     = (*x & mask_down<T>::value) * (*y & mask_down<T>::value);
        *(r+1) = (*x & mask_down<T>::value) * ((*y & mask_up<T>::value) >> (data_bits<T>::value/2));
    }

    template <typename T>
    inline void kernels_multiplication_block_up(T const* x, T const*  y, T * r)	
    {
        *r     = ((*x & mask_up<T>::value) >> (data_bits<T>::value/2) ) * (*y & mask_down<T>::value);	
        *(r+1) = ((*x & mask_up<T>::value) >> (data_bits<T>::value/2) ) * ((*y & mask_up<T>::value) >> (data_bits<T>::value/2));
    }

    template <typename T>
    inline void kernels_multiplication_base_reshaping(T const * a, T  const *  b, T * r)	
    {
        T q1 = (*(a+1) + *b) >> (data_bits<T>::value/2);
        T r1 = (*(a+1) + *b) & mask_down<T>::value;
        r1 *= base_half<T>::value;
        T q2 = (r1 + *a) >> data_bits<T>::value; 
        T r2 = (r1 + *a) & data_mask<T>::value;
        *r = r2;
        *(r+1) = q1 + q2 + *(b+1);
    }
    /**
    * END MULTIPLICATION BETWEEN TWO VLI
    */

    /**
    *
    * BEGINING MULTIPLICATION ALGO BETWEEN A VLI AND AN UNSIGNED INTEGER
    * 
    * This multiplication is done from right to left
    */
    template <typename BaseInt, std::size_t Size>
    void kernels_multiplication_number(BaseInt* x, BaseInt a)
    {
        BaseInt r[2] = {0,0};
        kernels_multiplication_block(&x[Size-1],&a,&(r[0]));
        x[Size-1] = r[0];
        for( std::size_t i = Size-1; i > 0; --i)
        {
            kernels_multiplication_block(&x[i-1],&a,&(r[0]));
            x[i-1] = r[0];
            x[i] += r[1];

            // Carry bit propagation
            for(std::size_t j = i; j < Size-2; ++j)
            {
                x[j+1] += x[j] >> data_bits<BaseInt>::value; //carry bit
                x[j] &= data_mask<BaseInt>::value; // Remove the carry bit
            }
        }
    }

    /**
    * END MULTIPLICATION BETWEEN TWO VLI
    */

   } //end namespace kernels
} //end namespace vli

#endif
