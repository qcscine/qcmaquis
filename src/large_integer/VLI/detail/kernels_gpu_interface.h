//
//  kernels_gpu_interface.h
//  vli
//
//  Created by Timothée Ewart on 22/08/11.
//  Copyright 2011 Université de Genève. All rights reserved.
//

#ifndef vli_kernels_gpu_interface_h
#define vli_kernels_gpu_interface_h

#include "kernels_gpu.h"
#include "detail/vli_size_param.hpp"

#include <cassert>

namespace vli {
namespace detail {


//NEW

/**
* global GPU functions declarations
*/
 
/**
* VLI_GPU functions
*/

template <typename T, int size>
__global__ void single_addition(T* x, T const* y);     

template <typename T>
__global__ void single_addition_int(T* x, T y);

template <typename T, int size>
__global__ void single_substraction(T* x, T const* y);  

template <typename T, int size>
__global__ void single_multiplication(T const* x, T const* y , T* z); 

template <typename T, int size>
__global__ void single_multiplication_int(T* x, int y);

template <typename T, int size>
__global__ void negate_gpu(T* x);


/**
* VLI_GPU_POLYNOMIAL functions
*/

template <typename T, int vli_size, int poly_size>
__global__ void monome_polynome_multiplication(T const* p, T const* m, T* res);

template <typename T, int vli_size, int poly_size>
__global__ void monome_polynome_addition(T* p, T const* m);

template <typename T, int vli_size, int poly_size> 
__global__ void polynome_polynome_addition(T* pa, T const* pb);

template <typename T, int vli_size, int poly_size>    
__global__ void polynome_polynome_substraction(T* pa, T const* pb); 

template <typename T, int vli_size, int poly_size>
__global__ void polynome_polynome_multiplication(T const* pa, T const* pb, T* pres);

/**
* VLI_GPU_VECTOR functions
*/

template <typename T, int vli_size, int poly_size>
__global__ void inner_prod_vector(T const* va, T const* vb, T* vres);

template <typename T, int vli_size, int poly_size>
__global__ void reduction_polynome(T const* va, T * vb, std::size_t SizeVector);

/**
* VLI_GPU functions
*/

void negate(type_vli::BaseInt* A) 
{ 
    dim3 dimgrid(1,1,1); 
    dim3 dimblock(1,1,1); 
    negate_gpu<type_vli::BaseInt, size_vli::value> <<< dimgrid, dimblock >>>(A); 
} 

void plus_assign(type_vli::BaseInt* A, type_vli::BaseInt const* B) 
{ 
    dim3 dimgrid(1,1,1);     
	dim3 dimblock(1,1,1); 
	single_addition<type_vli::BaseInt, size_vli::value> <<< dimgrid, dimblock >>>(A, B); 
} 

void plus_assign(type_vli::BaseInt* A,  type_vli::BaseInt a) 
{ 
    dim3 dimgrid(1,1,1);     
	dim3 dimblock(1,1,1); 
	single_addition_int<type_vli::BaseInt> <<< dimgrid, dimblock >>>(A, a); 
} 

void minus_assign(type_vli::BaseInt* A, type_vli::BaseInt const* B) 
{ 
    dim3 dimgrid(1,1,1); 
    dim3 dimblock(1,1,1); 
    single_substraction<type_vli::BaseInt, size_vli::value> <<< dimgrid, dimblock >>>(A, B); 
}     

void entrywise_multiplies(type_vli::BaseInt const* A, type_vli::BaseInt const* B, type_vli::BaseInt* C) 
{ 
    dim3 dimgrid(1,1,1); 
	dim3 dimblock(1,1,1); 
	single_multiplication<type_vli::BaseInt, size_vli::value> <<< dimgrid, dimblock >>>(A, B, C); 
} 

void multiplies_assign(type_vli::BaseInt* A, int a)
{
    dim3 dimgrid(1,1,1); 
	dim3 dimblock(1,1,1); 
    single_multiplication_int<type_vli::BaseInt,  size_vli::value> <<< dimgrid, dimblock >>>(A, a); 
}


/**
* VLI_GPU_POLYNOMIAL functions
*/
void poly_mono_multiply(type_vli::BaseInt const* A, type_vli::BaseInt const* B, type_vli::BaseInt* C)
{ 
    size_t size_poly_vli_value_square = size_poly_vli::value*size_poly_vli::value;
    dim3 dimgrid(1,1,1); 
    dim3 dimblock(size_poly_vli_value_square,1,1); 
    monome_polynome_multiplication<type_vli::BaseInt, size_vli::value, size_poly_vli::value>  <<< dimgrid, dimblock >>>(A, B, C); 
} 

void plus_assign_poly(type_vli::BaseInt* A, type_vli::BaseInt const* B) 
{ 
    size_t size_poly_vli_value_square = size_poly_vli::value*size_poly_vli::value;
    dim3 dimgrid(1,1,1);     
	dim3 dimblock(size_poly_vli_value_square,1,1); 
	polynome_polynome_addition<type_vli::BaseInt, size_vli::value, size_poly_vli::value> <<< dimgrid, dimblock >>>(A, B); 
} 

void minus_assign_poly(type_vli::BaseInt* A, type_vli::BaseInt const* B) 
{ 
    size_t size_poly_vli_value_square = size_poly_vli::value*size_poly_vli::value;
    dim3 dimgrid(1,1,1);     
	dim3 dimblock(size_poly_vli_value_square,1,1); 
	polynome_polynome_substraction<type_vli::BaseInt, size_vli::value, size_poly_vli::value> <<< dimgrid, dimblock >>>(A, B); 
} 

void poly_poly_multiply(type_vli::BaseInt const* A, type_vli::BaseInt const* B, type_vli::BaseInt* C) 
{ 
   	dim3 dimgrid(1,1,1); 
	dim3 dimblock(1,1,1); 
    polynome_polynome_multiplication<type_vli::BaseInt, size_vli::value, size_poly_vli::value> <<< dimgrid, dimblock >>>(A, B, C); 
} 

/**
* VLI_GPU_VECTOR functions
*/

void inner_product_vector(type_vli::BaseInt const* A, type_vli::BaseInt const* B, type_vli::BaseInt* C, std::size_t SizeVector, std::size_t threadsPerBlock) 
{ 
    std::size_t blocksPerGrid =  SizeVector/threadsPerBlock;
    inner_prod_vector<type_vli::BaseInt, size_vli::value, size_poly_vli::value> <<< blocksPerGrid,threadsPerBlock  >>>(A, B, C);  
} 
void vector_reduction(type_vli::BaseInt const* A, type_vli::BaseInt * B, std::size_t SizeVector) 
{ 
    //the reduction should be // if executed on one smp
    dim3 dimgrid(1,1,1); 
    dim3 dimblock(1,1,1); 
    reduction_polynome<type_vli::BaseInt, size_vli::value, size_poly_vli::value> <<< dimgrid, dimblock >>>(A, B, SizeVector);
}


/*
#define VLI_IMPLEMENT_GPU_KERNELS_FOR(r, data, TYPE) \
 \
void inner_prod_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size) \
{ \
    assert(false); \
} \

BOOST_PP_SEQ_FOR_EACH(VLI_IMPLEMENT_GPU_KERNELS_FOR, _, VLI_GPU_BASE_INT_TYPES_SEQ)

#undef VLI_IMPLEMENT_GPU_KERNELS_FOR
*/

}//detail
}//vli

#endif
