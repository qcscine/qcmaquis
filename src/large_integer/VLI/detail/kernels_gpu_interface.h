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
#include "detail/kernels_cpu_gpu.hpp" 
#include <cassert>

namespace vli {
namespace detail {

template <typename T>
__global__ void negate(T* x, int vli_size);

template <typename T>
__global__ void single_addition(T* x,   T const* y, int vli_size);     

template <typename T>
__global__ void single_int_addition(T* x, T y, int vli_size);

template <typename T>
__global__ void single_substraction(T* x,   T* y, int vli_size);  

template <typename T>    
__global__ void polynome_polynome_addition(T* x, T const* y,  int vli_size, int max_order );

template <typename T>    
__global__ void polynome_polynome_substraction(T* x, T const* y,  int vli_size, int max_order ); 

template <typename T>
__global__ void single_multiplication(T const* x, T const* y , T* z, int vli_size); 

template <typename T>
__global__ void single_multiplication_number( T *x, T a, int vli_size);

template <typename T>
__global__ void monome_polynome_multiplication(T const* p, T const* m, T* res, int vli_size, int max_order);

template <typename T>
__global__ void polynome_multication(T const* p1, T const* p2, T* res, int vli_size, int max_order);

template <typename T>
__global__ void inner_prod_vector(T const* p1, T const* p2, T* inter, int vli_size, int max_order, int size_vector);

template <typename T>
__global__ void reduction_polynome(T const* A, T * B,  int vli_size, int max_order, int size_vector);


template <class T, int size>
__device__ void kernels_addition_classic_gpu(T* x,  T const*  y)
{
//   using vli::detail::kernels_addition_classic;
//   kernels_addition_classic<T, size>(x,y);
}

template<typename T, int size>
__global__ void  add(T* A, T const* B)
{
     kernels_addition_classic_gpu<T, size>(A,B);
}

void plus_assign_lib(type_vli::BaseInt*  A, type_vli::BaseInt const*  B) 
{ 
    dim3 dimgrid(1,1,1);     
	dim3 dimblock(1,1,1); 
	add<type_vli::BaseInt, size_vli::value> <<< dimgrid, dimblock >>>(A, B); 
} 






#define VLI_IMPLEMENT_GPU_KERNELS_FOR(r, data, TYPE) \
void negate_gpu(TYPE* A, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    negate <<< dimgrid, dimblock >>>(A, vli_size); \
} \
void plus_assign_gpu(TYPE*  A, TYPE const*  B, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
	dim3 dimblock(1,1,1); \
	single_addition <<< dimgrid, dimblock >>>(A, B, vli_size); \
} \
void plus_assign_int_gpu(TYPE*  A, TYPE  B, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
	dim3 dimblock(1,1,1); \
	single_int_addition <<< dimgrid, dimblock >>>(A, B, vli_size); \
} \
void minus_assign_gpu(TYPE*  A, TYPE*  B, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    single_substraction <<< dimgrid, dimblock >>>(A, B, vli_size); \
}     \
void multiplies_assign_single_gpu(TYPE*  A, TYPE a, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    single_multiplication_number <<< dimgrid, dimblock >>>(A, a, vli_size); \
}     \
void entrywise_multiplies_gpu(TYPE const* a, TYPE const* b, TYPE* c, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
	dim3 dimblock(1,1,1); \
	single_multiplication <<< dimgrid, dimblock >>>(a, b , c, vli_size); \
} \
void inner_prod_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size) \
{ \
    assert(false); \
} \
void poly_multiply_gpu(TYPE const* a, TYPE const* b, TYPE* c, int vli_size, int max_order) \
{ \
   	dim3 dimgrid(1,1,1); \
	dim3 dimblock(1,1,1); \
    polynome_multication  <<< dimgrid, dimblock >>>(a, b , c, vli_size, max_order); \
} \
void poly_addition_gpu(TYPE* a, TYPE const* b, int vli_size, int max_order) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    polynome_polynome_addition  <<< dimgrid, dimblock >>>(a, b, vli_size, max_order); \
} \
void poly_substraction_gpu(TYPE* a, TYPE const* b, int vli_size, int max_order) \
{  \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    polynome_polynome_substraction  <<< dimgrid, dimblock >>>(a, b, vli_size, max_order); \
}     \
void poly_mono_multiply_gpu(TYPE const* a, TYPE const*b, TYPE* c, int vli_size, int max_order) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    monome_polynome_multiplication  <<< dimgrid, dimblock >>>(a, b, c ,vli_size, max_order); \
} \
void inner_product_vector_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size, int max_order, TYPE vector_size) \
{ \
    int threadsPerBlock = 16; \
    int blocksPerGrid = vector_size/16; \
    inner_prod_vector  <<< blocksPerGrid,threadsPerBlock  >>>(A, B, C,vli_size, max_order,vector_size);  \
} \
void vector_reduction_gpu(TYPE const* A, TYPE * B,  int vli_size, int max_order, TYPE vector_size) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    reduction_polynome <<< dimgrid, dimblock >>>(A, B, vli_size, max_order, vector_size); \
}

BOOST_PP_SEQ_FOR_EACH(VLI_IMPLEMENT_GPU_KERNELS_FOR, _, VLI_GPU_BASE_INT_TYPES_SEQ)

#undef VLI_IMPLEMENT_GPU_KERNELS_FOR


}//detail
}//vli

#endif
