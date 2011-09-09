//
//  kernels_gpu_interface.h
//  vli
//
//  Created by Timothée Ewart on 22/08/11.
//  Copyright 2011 Université de Genève. All rights reserved.
//

#include "detail/kernels_cpu_gpu.hpp"
#include "detail/vli_size_param.hpp"
#include <cassert>

//TO DO REMOVE ALL THESE REDICULOUS STATIC CAST INT -> SIZE_T


namespace vli {
namespace detail {

/**
* local declaration of the commun kernels
*/

template <typename BaseInt, std::size_t Size>
__host__ __device__ void kernels_addition_classic(BaseInt* x, BaseInt const* y);

template <typename BaseInt, std::size_t Size>
__host__ __device__ void kernel_negate_device(BaseInt* x);
 
template <typename BaseInt>
__host__ __device__ void kernels_addition_block(BaseInt* x, BaseInt const* y); 

template <typename BaseInt, std::size_t Size>
__host__ __device__ void kernels_addition_int(BaseInt* x, int a);

template <typename BaseInt, std::size_t Size>
__host__ __device__ void kernels_multiplication_classic_truncate(BaseInt * res, BaseInt const* x, BaseInt const* y);	

template <typename BaseInt>
__host__ __device__ void kernels_multiplication_block(BaseInt const* x, BaseInt const* y, BaseInt* r);

template <typename BaseInt>
__host__ __device__ void kernels_multiplication_block_down(BaseInt const* x, BaseInt const*  y, BaseInt* r);

template <typename BaseInt>
__host__ __device__ void kernels_multiplication_block_up(BaseInt const* x, BaseInt const*  y, BaseInt* r);	

template <typename BaseInt>
__host__ __device__ void kernels_multiplication_base_reshaping(BaseInt const* x, BaseInt  const*  y, BaseInt* r);	

template <typename BaseInt, std::size_t Size>
__host__ __device__ void kernels_multiplication_number(BaseInt* x, BaseInt a);

template <typename BaseInt, int size>
__device__ void single_multiplication_device(BaseInt const* x, BaseInt const* y, BaseInt* z);  

template <typename T, int vli_size, int poly_size>
__device__ void polynome_polynome_multiplication_device(T const* p1, T const* p2, T* res);


/**
* a kind of hook function with a little bit of arithmetic in case of signed int (multiplication)
*/

/**
* VLI_GPU functions
*/

template <typename T, int size>
__global__ void negate_gpu(T* x)
{
     kernel_negate_device<T,static_cast<std::size_t>(size)>(x);
}

template <typename T, int size>
__global__ void single_addition(T* x, T const* y)     
{
    kernels_addition_classic<T,static_cast<std::size_t>(size)>(x,y);
}    

template <typename T>
__global__ void single_addition_int(T* x, T y)
{
    //TODO BUG: what about the carry bit?
    kernels_addition_block<T>(x,&y);  
}

template <typename T, int size>
__global__ void single_substraction(T* x, T const* y)
{
    kernel_negate_device<T,static_cast<std::size_t>(size)>(const_cast<T*>(y));
    kernels_addition_classic<T,static_cast<std::size_t>(size)>(x,y);
}

//__global__ and __device__ is impossible in the same time
template <typename T, int size>
__global__  void single_multiplication(T const* x, T const* y, T* z)     
{
     single_multiplication_device<T,size>(x,y,z);
}

template <typename T, int size>
__device__  void single_multiplication_device(T const* x, T const* y, T* z)  
{
    int na(1),nb(1);
/*
    bool result_is_negative = static_cast<bool>((x[size-1] ^ y[size-1]) >> data_bits<T>::value);
    if(result_is_negative)// test if 
    {
        kernel_negate_device<T,static_cast<std::size_t>(size)>(const_cast<T* >(x));
        kernels_multiplication_classic_truncate<T,static_cast<std::size_t>(size)>(z,x,y); 
        kernel_negate_device<T,static_cast<std::size_t>(size)>(const_cast<T* >(x));
    }
    else
    {
        kernels_multiplication_classic_truncate<T,static_cast<std::size_t>(size)>(z,x,y); 
    }
*/
       
    if( static_cast<bool>((x[size-1]) >> data_bits<T>::value)){
        kernel_negate_device<T,static_cast<std::size_t>(size)>(const_cast<T* >(x));
        na = -1;
    }

    if( static_cast<bool>((y[size-1]) >> data_bits<T>::value)){
        kernel_negate_device<T,static_cast<std::size_t>(size)>(const_cast<T* >(y));
        nb = -1;
    }

    kernels_multiplication_classic_truncate<T,static_cast<std::size_t>(size)>(z,x,y);

    if(nb*na == -1)
        kernel_negate_device<T,static_cast<std::size_t>(size)>(z);
       
    if(na == -1){
        kernel_negate_device<T,static_cast<std::size_t>(size)>(const_cast<T* >(x));
    }

    if(nb == -1){
        kernel_negate_device<T,static_cast<std::size_t>(size)>(const_cast<T* >(y));
    }
}

template <typename T, int size>
__global__ void single_multiplication_int(T* x, int y)
{
    int na(1),nb(1);
    
    if( static_cast<bool>((x[size-1]) >> data_bits<T>::value)){
        kernel_negate_device<T,static_cast<std::size_t>(size)>(x);
        na = -1;
    }
    
    if(y<0){
        y *=-1;
        nb = -1;
    }
     
    kernels_multiplication_number<T,static_cast<std::size_t>(size)>(x,y); 
     
    if(nb*na == -1)
        kernel_negate_device<T,static_cast<std::size_t>(size)>(x);
}

/**
I try to implement this stuff into the commun kernel impossible ... (3 days of trying)
the compiler does not want, we should move it inside ! 
**/
template <typename BaseInt, std::size_t Size>
__device__ void kernel_negate_device(BaseInt* x)
{
    BaseInt one(1);
    #pragma unroll
    for(int i=0; i < Size-1; ++i)
        *(x+i) = (~*(x+i))&data_mask<BaseInt>::value;
    *(x+Size-1) = (~*(x+Size-1))&(base<BaseInt>::value+data_mask<BaseInt>::value);
    
    kernels_addition_block(x,&one);
}

/**
* VLI_GPU_MONOMIAL functions
*/

template <typename T, int vli_size, int poly_size> 
__global__ void monome_polynome_multiplication(T const* p, T const* m, T* res, std::size_t j_exp, std::size_t h_exp)
{ 
    std::size_t max_order = static_cast<std::size_t>(poly_size);
    std::size_t offset0(0);
    std::size_t offset1(0);
    
    for(std::size_t je = 0; je < max_order-j_exp; ++je)
    {
        for(std::size_t he = 0; he < max_order-h_exp; ++he)
        {   
            offset0 = (je*max_order+he)*static_cast<std::size_t>(vli_size);
            offset1 = ((j_exp+je)*max_order+h_exp+he)*static_cast<std::size_t>(vli_size);
            single_multiplication_device<T,static_cast<std::size_t>(vli_size)>(&p[offset0],&m[0],&res[offset1]);        
        }
    }

}  

/**
* VLI_GPU_POLYNOMIAL functions
*/

template <typename T, int vli_size, int poly_size> 
__global__ void polynome_polynome_addition(T* x, T const* y) 
{   
    std::size_t xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
	std::size_t offset = xIndex*vli_size;    
    kernels_addition_classic<T,static_cast<std::size_t>(vli_size)>(&x[offset],&y[offset]);    //1 see line 148
}
  

template <typename T, int vli_size, int poly_size> 
__global__ void polynome_polynome_substraction(T* x, T const* y) 
{
    std::size_t xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
	std::size_t offset = xIndex*vli_size;    
    kernel_negate_device<T,static_cast<std::size_t>(vli_size)>(const_cast<T*>(&y[offset]));
    kernels_addition_classic<T,static_cast<std::size_t>(vli_size)>(&x[offset],&y[offset]);
}


template <typename T, int vli_size, int poly_size>
__global__ void polynome_polynome_multiplication(T const* x, T const* y, T* res)
{
    polynome_polynome_multiplication_device<T, vli_size, poly_size>(x,y,res);
}

template <typename T, int vli_size, int poly_size>
__device__ void polynome_polynome_multiplication_device(T const* p1, T const* p2, T* res)
{
    for(std::size_t je1 = 0; je1 < poly_size; ++je1)
    {
        for(std::size_t he1 = 0; he1 < poly_size; ++he1)
        {
            for(std::size_t je2 = 0; je2 < poly_size - je1; ++je2)
            {
                for(std::size_t he2 = 0; he2 < poly_size - he1; ++he2)
                {
                    T inter[vli_size];
                    #pragma unroll
                    for(std::size_t i=0 ; i < vli_size ;++i)
                        inter[i] = 0;
                    //memset(&inter[0], 0 , vli_size*sizeof(T)); memset forbiden on GPU it is corrupted the memory, not stable !!!
                    std::size_t offset0 = ((je1+je2)*poly_size + he1+he2)*vli_size;
                    std::size_t offset1 = (je1*poly_size+he1)*vli_size;                    
                    std::size_t offset2 = (je2*poly_size+he2)*vli_size;
                    single_multiplication_device<T,static_cast<std::size_t>(vli_size)>(&p1[offset1],&p2[offset2],&inter[0]);        
                    kernels_addition_classic<T,static_cast<std::size_t>(vli_size)>(&res[offset0],&inter[0]);
                } 
            }
        }      
    }
} 

/**
* VLI_GPU_VECTOR functions
*/

    
template  <typename T, int vli_size, int poly_size>
__global__ void inner_prod_vector(T const* v1, T const* v2, T* res)
{
    unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
    unsigned int size_poly = vli_size*poly_size*poly_size;
    unsigned int offset = xIndex*size_poly;
    //multiplication between polynomial
    polynome_polynome_multiplication_device<T,vli_size,poly_size>(&v1[offset],&v2[offset],&res[offset]); 
}
    
    
template  <typename T, int vli_size, int poly_size>
__global__ void reduction_polynome(T const* v1, T* v2, std::size_t SizeVector)
{ 
    /*
    *  WARNING the kernels can only be instancied with one value
    *  here vli_size, so, one more loop.
    *  SOLUTION : BOOST_PP_SEQ_FOR_EACH ..... 
    */    
    std::size_t size_poly = static_cast<std::size_t>(vli_size*poly_size*poly_size);  
    std::size_t offset0(0);  
    std::size_t offset1(0);  
    
    for(std::size_t i=0 ; i < SizeVector ; ++i){
        for(std::size_t j=0 ; j < poly_size*poly_size ; ++j){ //additional loop
            offset0 = j*vli_size;
            offset1 = i*size_poly+j*vli_size;
            kernels_addition_classic<T,vli_size>(&v2[offset0],&v1[offset1]);
        }
    }
    /*
     for(std::size_t i=0 ; i < SizeVector ; ++i){
        kernels_addition_classic<T,vli_size*poly_size*poly_size>(&v2[i],&v1[i]); 
         
     }*/
    
    
}
   
template <typename T, int vli_size>
__global__ void polynome_int_addition(T* x, int y)
{
    kernels_addition_int<T, vli_size>(x,y);
}

/**
  * The C++ functions that call the kernels
  */
template <typename BaseInt, std::size_t Size>
void negate(BaseInt* A)
{
    dim3 dimgrid(1,1,1);
    dim3 dimblock(1,1,1);
    negate_gpu<BaseInt, Size> <<< dimgrid, dimblock >>>(A);
}

template <typename BaseInt, std::size_t Size>
void plus_assign(BaseInt* A, BaseInt const* B)
{
    dim3 dimgrid(1,1,1);
	dim3 dimblock(1,1,1);
	single_addition<BaseInt, Size> <<< dimgrid, dimblock >>>(A, B);
}

template <typename BaseInt, std::size_t Size>
void plus_assign(BaseInt* A, BaseInt a)
{
    dim3 dimgrid(1,1,1);
	dim3 dimblock(1,1,1);
	single_addition_int<BaseInt> <<< dimgrid, dimblock >>>(A, a);
}

template <typename BaseInt, std::size_t Size>
void minus_assign(BaseInt* A, BaseInt const* B)
{
    dim3 dimgrid(1,1,1);
    dim3 dimblock(1,1,1);
    single_substraction<BaseInt, Size> <<< dimgrid, dimblock >>>(A, B);
}

template <typename BaseInt, std::size_t Size>
void entrywise_multiplies(BaseInt const* A, BaseInt const* B, BaseInt* C)
{
    dim3 dimgrid(1,1,1);
	dim3 dimblock(1,1,1);
	single_multiplication<BaseInt, Size> <<< dimgrid, dimblock >>>(A, B, C);
} 

template <typename BaseInt, std::size_t Size>
void multiplies_assign(BaseInt* A, int a)
{
    //TODO is 'int a' ok?
    dim3 dimgrid(1,1,1);
	dim3 dimblock(1,1,1);
    single_multiplication_int<BaseInt, Size> <<< dimgrid, dimblock >>>(A, a);
}


/**
* VLI_GPU_POLYNOMIAL functions
*/
template <typename BaseInt, std::size_t Size>
void poly_mono_multiply(BaseInt const* A, BaseInt const* B, BaseInt* C, std::size_t j_exp, std::size_t h_exp)
{
    //size_t size_poly_vli_value_square = size_poly_vli::value*size_poly_vli::value;
    dim3 dimgrid(1,1,1);
    dim3 dimblock(1,1,1);
    // TODO size_poly
    monome_polynome_multiplication<BaseInt, Size, size_poly_vli::value>  <<< dimgrid, dimblock >>>(A, B, C, j_exp, h_exp);
}

template <typename BaseInt, std::size_t Size>
void plus_assign_poly_int(BaseInt* A, int a)
{
    // TODO is 'int a' ok?
    dim3 dimgrid(1,1,1);
    dim3 dimblock(1,1,1);
    polynome_int_addition<BaseInt, Size> <<< dimgrid, dimblock>>>(A,a);
}

template <typename BaseInt, std::size_t Size>
void plus_assign_poly(BaseInt* A, BaseInt const* B)
{
    dim3 dimgrid(1,1,1);
	dim3 dimblock(size_poly_vli::value*size_poly_vli::value,1,1);
    // TODO size_poly 
	polynome_polynome_addition<BaseInt, Size, size_poly_vli::value> <<< dimgrid, dimblock >>>(A, B);
}

template <typename BaseInt, std::size_t Size>
void minus_assign_poly(BaseInt* A, BaseInt const* B)
{
    size_t size_poly_vli_value_square = size_poly_vli::value*size_poly_vli::value;
    dim3 dimgrid(1,1,1);
	dim3 dimblock(size_poly_vli_value_square,1,1);
    //TODO size_poly
	polynome_polynome_substraction<BaseInt, Size, size_poly_vli::value> <<< dimgrid, dimblock >>>(A, B);
}

template <typename BaseInt, std::size_t Size>
void poly_poly_multiply(BaseInt const* A, BaseInt const* B, BaseInt* C)
{
  	dim3 dimgrid(1,1,1);
    dim3 dimblock(1,1,1); 
    //TODO size_poly
    polynome_polynome_multiplication<BaseInt, Size, size_poly_vli::value> <<< dimgrid, dimblock >>>(A, B, C); 
} 

/**
* VLI_GPU_VECTOR functions
*/

template <typename BaseInt, std::size_t Size>
void inner_product_vector(BaseInt const* A, BaseInt const* B, BaseInt* C, std::size_t vector_size, std::size_t threads_per_block) 
{
    std::size_t blocksPerGrid =  vector_size/threads_per_block;
    //TODO size_poly
    inner_prod_vector<BaseInt, Size, size_poly_vli::value> <<< blocksPerGrid,threads_per_block  >>>(A, B, C);  
} 

template <typename BaseInt, std::size_t Size>
void vector_reduction(BaseInt const* A, BaseInt * B, std::size_t vector_size)
{
    //the reduction should be // if executed on one smp
    dim3 dimgrid(1,1,1);
    dim3 dimblock(1,1,1);
    reduction_polynome<BaseInt, Size, size_poly_vli::value> <<< dimgrid, dimblock >>>(A, B, vector_size);
}

} // namespace detail
} // namespace vli

#include "kernels_gpu.h"

