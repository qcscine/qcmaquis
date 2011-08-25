//
//  kernels_gpu_interface.h
//  vli
//
//  Created by Timothée Ewart on 22/08/11.
//  Copyright 2011 Université de Genève. All rights reserved.
//

#ifdef VLI_GPU_DEBUG
#include <cstdio>
#endif //VLI_GPU_DEBUG
//#include "kernels_gpu.h"
#include "detail/kernels_cpu_gpu.hpp" 
#include "kernels_gpu_interface.h"
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
__host__ __device__ void kernels_multiplication_classic_truncate(BaseInt * res, BaseInt const* x, BaseInt const* y);	

template <typename BaseInt>
__host__ __device__  void kernels_multiplication_block(BaseInt const* x, BaseInt const* y, BaseInt* r);

template <typename BaseInt>
__host__ __device__ void kernels_multiplication_block_down(BaseInt const* x, BaseInt const*  y, BaseInt* r);

template <typename BaseInt>
__host__ __device__ void kernels_multiplication_block_up(BaseInt const* x, BaseInt const*  y, BaseInt* r);	

template <typename BaseInt>
__host__ __device__ void kernels_multiplication_base_reshaping(BaseInt const* x, BaseInt  const*  y, BaseInt* r);	

template <typename BaseInt, std::size_t Size>
__host__ __device__ void kernels_multiplication_number(BaseInt* x, BaseInt a);

template <typename BaseInt, int size>
__device__  void single_multiplication_device(BaseInt const* x, BaseInt const* y, BaseInt* z);  

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
__global__ void monome_polynome_multiplication(T const* p, T const* m, T* res)
{
    std::size_t xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
	std::size_t offset = xIndex*vli_size;    
    single_multiplication_device<T,static_cast<std::size_t>(vli_size)>(&p[offset],&m[0],&res[offset]);        
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
                    memset(&inter[0], 0 , vli_size*sizeof(T));
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

  
/*    
     
    
    
template <typename T>
__global__ void inner_prod_vector(T const* p1, T const* p2, T* inter, int vli_size, int max_order, int size_vector)
{
    unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
    unsigned int size_poly = vli_size*max_order*max_order;
    unsigned int offset = xIndex*size_poly;
        //mutiplication
    polynome_polynome_multiplication(&p1[offset],&p2[offset],&inter[offset],vli_size,max_order); 
}
    
    
template <typename T>
__global__ void reduction_polynome(T const* A, T * B,  int vli_size, int max_order, int size_vector)
{ 
    int size_poly = vli_size*max_order*max_order;
    
    for(unsigned int i=0 ; i < size_vector ; i++)
            addition_classic_kernel_gpu(&B[0],&A[i*size_poly], size_poly);  
            
}
    
    
template <typename T>
__device__ void copy_kernel_gpu(T* x, T const* y)
{
    *x = *y;
}    
    

template <typename T>
__device__ void addition_with_int_kernel_gpu(T* x, int y, int vli_size)
{
    *x  += y;

    for (int k = 0; k < vli_size-1; ++k)
    { 
        *(x+k+1)  += *(x+k) >> data_bits<T>::value;
        *(x+k)    &= data_mask<T>::value;
    }
    *(x+vli_size-1) &= base<T>::value + data_mask<T>::value;
}



template <typename T>
__device__ void addition_classic_kernel_gpu(T* x, T const* y, int vli_size)
{

    for (int k = 0; k < vli_size-1; ++k)
    { 
        *(x+k)    += *(y+k);
        *(x+k+1)  += *(x+k) >> data_bits<T>::value;
        *(x+k)    &= data_mask<T>::value;
    }
    
    *(x+vli_size-1) += *(y+vli_size-1);
    *(x+vli_size-1) &= base<T>::value + data_mask<T>::value;
}    
    


*/
}//detail
}//vli


