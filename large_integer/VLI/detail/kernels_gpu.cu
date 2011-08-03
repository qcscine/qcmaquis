
/**
    Do not forget NVCC will compile this file, it dislikes STL and boost ....
    Just very basic C++ or C (template ok). 
 
	The link between the DMRG code (CPU) and the DMRG code (GPU) is done inside the file  kernels_gpu.h
	where the definition of function/wrapper are presented
*/

#ifdef VLI_GPU_DEBUG
#include <cstdio>
#endif //VLI_GPU_DEBUG

#include "kernels_gpu.h"
#include "detail/bit_masks.hpp"
#include <cassert>


namespace vli {
namespace detail {

template <typename T>
__device__ void copy_kernel_gpu(T* x, T const* y)
{
    *x = *y;
}    
    
template <typename T>
__device__ void negate_device(T* x, int vli_size)
{
    #pragma unroll
    for(int i=0; i < vli_size-1; ++i)
        *(x+i) = (~*(x+i))&data_mask<T>::value;

    *(x+vli_size-1) = (~*(x+vli_size-1))&(base<T>::value+data_mask<T>::value);
    
    addition_with_int_kernel_gpu(x, 1, vli_size);
}

template <typename T>
__device__ void addition_with_int_kernel_gpu(T* x, int y, int vli_size)
{
    *x    += y;
  
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

template <typename T>
__global__ void negate(T* x, int vli_size)
{
    negate_device(x,vli_size);
}

template <typename T>
__global__ void single_addition(T* x,   T const* y, int vli_size)     
{
    addition_classic_kernel_gpu(x, y, vli_size);    
}    
    

template <typename T>
__global__ void single_substraction(T* x,   T* y, int vli_size)     
{
    negate_device(y,vli_size);
    addition_classic_kernel_gpu(x, y, vli_size);        
}    

    

template <typename T>    
__global__ void polynome_polynome_addition(T* x, T const* y,  int vli_size, int max_order ) 
{
    const int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
	int offset(0);
    
	if(xIndex < 2) //useless, because we set the number of thread to one
	{
        for(int i=0; i< max_order*max_order;++i){
            addition_classic_kernel_gpu(&x[offset],&y[offset],vli_size);    //1 see line 148
            offset += vli_size;
        }
    }
}
    
template <typename T>    
__global__ void polynome_polynome_substraction(T* x, T const* y,  int vli_size, int max_order ) 
{
    const int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
    int offset(0);
        
    if(xIndex < 2) //useless, because we set the number of thread to one
    {
        for(int i=0; i< max_order*max_order;++i){
            negate_device(&x[offset],vli_size);
            addition_classic_kernel_gpu(&x[offset],&y[offset],vli_size);    //1 see line 148
            offset += vli_size;
        }
    }
}
        
    
/**
classical multiplication, the operation " addition " is serial but we sum all number together 
*/
template <typename T>
__device__ void addition_kernel_gpu_noiter(T* x, const T* y )
{
	*(x)    +=  *(y);
	*(x+1)  += *(x) >> data_bits<T>::value;
	*(x)    &= data_mask<T>::value;
}

template <typename T>
__device__ void multiplication_kernel_up_gpu(const T*   x, const T*    y, T * r)	
{
	*r	    = ((*x & mask_up<T>::value) >> (data_bits<T>::value/2) ) * (*y & mask_down<T>::value);	
	*(r+1)	= ((*x & mask_up<T>::value) >> (data_bits<T>::value/2) ) * ((*y & mask_up<T>::value) >> (data_bits<T>::value/2));
}		
	
template <typename T>
__device__ void multiplication_kernel_down_gpu(const T*  x,  const T*   y, T * r)	
{	
	*r     = (*x & mask_down<T>::value) * (*y & mask_down<T>::value);
	*(r+1) = (*x & mask_down<T>::value) * ((*y & mask_up<T>::value) >> (data_bits<T>::value/2));
}

template <typename T>
__device__ void multiplication_kernel_base_reshaping_gpu(T*  a, const T*   b, T * r)	
{	
	T q1 = (*(a+1) + *b) >> (data_bits<T>::value/2);
	T r1 = (*(a+1) + *b) & mask_down<T>::value;
	r1 *= base_half<T>::value;
	T q2 = (r1 + *a) >> data_bits<T>::value;
	T r2 = (r1 + *a) & data_mask<T>::value;
	*r  = r2;
	*(r+1) = q1 + q2 + *(b+1);
}

template <typename T>
__device__  void multiplication_block_gpu( const T*   x,   const  T*   y, T *r)	
{
    T a[2] = {0,0};
	T b[2] = {0,0};
	/**
	 Divide and conquer algo (see my notes - for the euclidian division tips)
		X <=> Xl Xr (half of the binary number)
	 x  Y <=> Yl Yr (half of the binary number)
	-------
	= 2^n XlYl + 2^(n/2) (XlYr + XrYl) + XrYr (multiplication_kernel_gpu_down and multiplication_kernel_gpu_up)
	------- 
	= (q1+q2 + Xl*Yl)*BASE + r2  (multiplication_kernel_base_reshaping)
	*/
	multiplication_kernel_down_gpu(x,y,a);
	multiplication_kernel_up_gpu(x,y,b);
	multiplication_kernel_base_reshaping_gpu(a,b, r);
}	

template <typename T>
__device__ void multiplication_classic_kernel_gpu(const T* x,  const T* y , T* z, int vli_size)// remove num_int maybe ?
{
	T r[2] = {0,0};	//for local block calculation
	
    for (int i = 0 ; i < vli_size; ++i) 
    {
        for(int j = 0 ; j < vli_size ; ++j)  			
        {	
            int m = j + i;
            multiplication_block_gpu(&x[i], &y[j], &r[0]);
            addition_kernel_gpu_noiter(&z[m] , &r[0]);
            addition_kernel_gpu_noiter(&z[m+1], &r[1]);
        }
    }
}
    
template <typename T>
__global__ void single_multiplication(T const* x, T const* y , T* z, int vli_size)     
{
//   bool result_is_negative = static_cast<bool>((x[Size-1] ^ y[Size-1]) >> data_bits<BaseInt>::value);
//     if(result_is_negative)// test if 
//   {
//      negate_device(y,vli_size);
        multiplication_classic_kernel_gpu(x, y, z, vli_size);  // +*- or -*+ 
//      negate_device(y,vli_size);
//   }
//   else
//   {
//       multiplication_classic_kernel_gpu(x, y, z, vli_size);  // +*+ or -*-  
//   }
}
    
template <typename T>
__global__ void monome_polynome_multiplication(T const* p, T const* m, T* res, int vli_size, int max_order)
{
    const int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
	
    #pragma unroll
    for(std::size_t i = 0 ; i < max_order*max_order ; i++)
    {
        std::size_t offset = i*vli_size;
        multiplication_classic_kernel_gpu(&p[offset],&m[0],&res[offset],vli_size);
    }
    
}
   
template <typename T>
__device__ void polynome_polynome_multiplication(T const* p1, T const* p2, T* res, int vli_size, int max_order)
{
    for(std::size_t je1 = 0; je1 < max_order; ++je1)
    {
        for(std::size_t he1 = 0; he1 < max_order; ++he1)
        {
            for(std::size_t je2 = 0; je2 < max_order - je1; ++je2)
            {
                for(std::size_t he2 = 0; he2 < max_order - he1; ++he2)
                {
                    T inter[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // to do find better
                    std::size_t offset0 = ((je1+je2)*max_order + he1+he2)*vli_size;
                    std::size_t offset1 = (je1*max_order+he1)*vli_size;                    
                    std::size_t offset2 = (je2*max_order+he2)*vli_size;
                    multiplication_classic_kernel_gpu(&p1[offset1],&p2[offset2],&inter[0],vli_size);
                    addition_classic_kernel_gpu(&res[offset0],&inter[0],vli_size);
                } 
            }
        }      
    }
} 

template <typename T>
__global__ void polynome_multication(T const* p1, T const* p2, T* res, int vli_size, int max_order)
{
    polynome_polynome_multiplication(p1,p2,res,vli_size,max_order);
}
    
template <typename T>
__global__ void inner_prod_vector(T const* p1, T const* p2, T* inter, int vli_size, int max_order, int size_vector)
{
    unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
    unsigned int size_poly = vli_size*max_order*max_order;
    unsigned int offset = xIndex*size_poly;
        //mutiplication
    polynome_polynome_multiplication(&p1[offset],&p2[offset],&inter[offset],vli_size,max_order); 
}
    
    
/**
   the reduction is serial due to the race conditions   
*/
template <typename T>
__global__ void reduction_polynome(T const* A, T * B,  int vli_size, int max_order, int size_vector)
{ 
    unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; 
    int size_poly = vli_size*max_order*max_order;
    
    for(unsigned int i=0 ; i < size_vector ; i++)
            addition_classic_kernel_gpu(&B[0],&A[i*size_poly], size_poly);   
}
    

#define VLI_IMPLEMENT_GPU_KERNELS_FOR(r, data, TYPE) \
void negate_gpu(TYPE* A, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    negate <<< dimgrid, dimblock >>>(A, vli_size); \
} \
void plus_assign_gpu(TYPE*  A, TYPE const*  B, int num_integers, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
	dim3 dimblock(1,1,1); \
	single_addition <<< dimgrid, dimblock >>>(A, B, vli_size); \
} \
void minus_assign_gpu(TYPE*  A, TYPE*  B, int vli_size) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    single_substraction <<< dimgrid, dimblock >>>(A, B, vli_size); \
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
void inner_product_vector_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size, int max_order, int vector_size) \
{ \
    int threadsPerBlock = 512; \
    int blocksPerGrid = vector_size/512; \
    inner_prod_vector  <<< blocksPerGrid,threadsPerBlock  >>>(A, B, C,vli_size, max_order,vector_size);  \
} \
void vector_reduction_gpu(TYPE const* A, TYPE * B,  int vli_size, int max_order, int vector_size) \
{ \
    dim3 dimgrid(1,1,1); \
    dim3 dimblock(1,1,1); \
    reduction_polynome <<< dimgrid, dimblock >>>(A, B, vli_size, max_order, vector_size); \
}

BOOST_PP_SEQ_FOR_EACH(VLI_IMPLEMENT_GPU_KERNELS_FOR, _, VLI_GPU_BASE_INT_TYPES_SEQ)

#undef VLI_IMPLEMENT_GPU_KERNELS_FOR

}//detail
}//vli


