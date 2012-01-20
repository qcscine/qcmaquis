//
//  kernels_gpu_interface.h
//  vli
//
//  Created by Timothée Ewart on 22/08/11.
//  Copyright 2011 Université de Genève. All rights reserved.
//

#include "vli/detail/kernels_cpu_gpu.hpp"
#include "vli/detail/kernels_gpu.h"


#include "utils/cuPrintf.cu"

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
    __host__ __device__ void kernels_multiplication_block(BaseInt const* x, BaseInt const* y, BaseInt* r);
 
    template <typename BaseInt>
    __host__ __device__ void kernels_multiplication_block_down(BaseInt const* x, BaseInt const*  y, BaseInt* r);

    template <typename BaseInt>
    __host__ __device__ void kernels_multiplication_block_up(BaseInt const* x, BaseInt const*  y, BaseInt* r);	

    template <typename BaseInt>
    __host__ __device__ void kernels_multiplication_base_reshaping(BaseInt const* x, BaseInt  const*  y, BaseInt* r);	

    template <typename BaseInt, std::size_t Size>
    __device__ void single_multiplication_device(BaseInt const* x, BaseInt const* y, BaseInt* z);  

    template <typename BaseInt, std::size_t Size>
    __device__ void polynome_polynome_multiplication_device(BaseInt const* p1, BaseInt const* p2, BaseInt* res);

// functions for the blocks algos

    template <typename BaseInt, std::size_t Size>
    __device__ void algo_triangle_up(int block_ai, int block_bj,unsigned int Order, BaseInt const* a,  BaseInt const* b, BaseInt* c);

    template <typename BaseInt, std::size_t Size>
    __device__ void algo_triangle_down(int block_ai, int block_bj,unsigned int Order, BaseInt const* a,  BaseInt const* b, BaseInt* c);
        
    template <typename BaseInt, std::size_t Size>
    __device__ void algo_diag(int block_ai, int block_bj,unsigned int Order, BaseInt const* a, BaseInt const* b, BaseInt* c);
    
    template <typename BaseInt, std::size_t Size>
    __device__ void algo_block_algo(int i, int j,unsigned int Order, BaseInt const* a, BaseInt const* b, BaseInt* c);
        
/**
* a kind of hook function with a little bit of arithmetic in case of signed int (multiplication)
*/

/**
* VLI_GPU functions
*/

template <typename BaseInt, std::size_t Size>
__device__  void single_multiplication_device(BaseInt const* x, BaseInt const* y, BaseInt* z)  
{
 //   int na(1),nb(1);

    bool result_is_negative = static_cast<bool>((x[Size-1] ^ y[Size-1]) >> data_bits<BaseInt>::value);
    if(result_is_negative)// test if 
    {
        kernel_negate_device<BaseInt,Size>(const_cast<BaseInt* >(x));
        kernels_multiplication_classic_truncate<BaseInt,Size>(z,x,y); 
        kernel_negate_device<BaseInt,Size>(const_cast<BaseInt* >(x));
    }
    else
    {
        kernels_multiplication_classic_truncate<BaseInt,Size>(z,x,y); 
    }
}
/**
I try to implement this stuff into the commun kernel impossible ... (3 days of trying)
the compiler does not want, we should move it inside ! 
**/
template <typename BaseInt, std::size_t Size>
__device__ void kernel_negate_device(BaseInt* x)
{
    BaseInt one(1);
    for(std::size_t i=0; i < Size-1; ++i)
        *(x+i) = (~*(x+i))&data_mask<BaseInt>::value;
    *(x+Size-1) = (~*(x+Size-1))&(base<BaseInt>::value+data_mask<BaseInt>::value);
    
    // TODO BUG!!!!
    kernels_addition_block(x,&one);
}
    
template <typename BaseInt, std::size_t Size>
__device__ void polynome_polynome_multiplication_device(unsigned int max_order, BaseInt const* p1, BaseInt const* p2, BaseInt* res)
{
    for(std::size_t je1 = 0; je1 < max_order; ++je1)
    {               
        for(std::size_t je2 = 0; je2 < max_order; ++je2)
        {
            for(std::size_t he1 = 0; he1 < max_order; ++he1)
            {                
                for(std::size_t he2 = 0; he2 < max_order; ++he2)
                {            
                    BaseInt inter[Size];
                    #pragma unroll
                    for(std::size_t i=0 ; i < Size ;++i)
                        inter[i] = 0;
        
                    std::size_t offset0 = ((je1+je2)*2*max_order + he1+he2)*Size;
                    std::size_t offset1 = (je1*max_order+he1)*Size;
                    std::size_t offset2 = (je2*max_order+he2)*Size;
                    single_multiplication_device<BaseInt,Size>(&p1[offset1],&p2[offset2],&inter[0]);
                    kernels_addition_classic<BaseInt,Size>(&res[offset0],&inter[0]);
                } 
            }
        }
    }
} 
        
        
template <typename BaseInt, std::size_t Size>
void algo_triangle_up(int block_ai, int block_bj, unsigned int Order, BaseInt const* a,  BaseInt const* b, BaseInt *c){
    std::size_t n(0);
    int offset_block_ai = block_ai*Order;
    int offset_block_bj = block_bj*Order;            
    int offset_a, offset_b, offset_c;
    
    for(int i = 0; i < Order-1; ++i){
        for(int j = n; j < Order-1; ++j){
            BaseInt inter[Size];
            #pragma unroll
            for(std::size_t k=0 ; k < Size ;++k)
                inter[k] = 0;
            
            offset_a = (j-n+offset_block_ai)*Size;
            offset_b = (n+offset_block_bj)*Size;
            offset_c = ((offset_block_ai+offset_block_bj)*2+j)*Size;
             
            single_multiplication_device<BaseInt,Size>(&a[offset_a],&b[offset_b],inter);
            kernels_addition_classic<BaseInt,Size>(&c[offset_c],inter);
        }    
        n++;
    }        
}

template <typename BaseInt, std::size_t Size>
void algo_diag(int block_ai, int block_bj,  unsigned int Order, BaseInt const* a, BaseInt const* b, BaseInt *c){
    int OrderMinusOne = Order-1;
    int offset_block_ai = block_ai*Order;
    int offset_block_bj = block_bj*Order;   
    int offset_a, offset_b, offset_c;
    
    for(int i = 0; i < Order; ++i){
        BaseInt inter[Size];
        #pragma unroll
        for(std::size_t k=0 ; k < Size ;++k)
            inter[k] = 0;
        
        offset_a = (offset_block_bj+Order-1-i)*Size;
        offset_b = (offset_block_ai+i)*Size;
        offset_c = ((offset_block_ai+offset_block_bj)*2+OrderMinusOne)*Size;

        single_multiplication_device<BaseInt,Size>(&a[offset_a],&b[offset_b],inter);
        kernels_addition_classic<BaseInt,Size>(&c[offset_c],inter);
    }
    
}

template <typename BaseInt, std::size_t Size>
void algo_triangle_down(int block_ai, int block_bj, unsigned int Order, BaseInt const* a, BaseInt const* b, BaseInt *c){
    int n(0);
    int offset_block_ai = (block_ai+1)*(Order)-1;
    int offset_block_bj = (block_bj+1)*(Order)-1;
    int offset_a, offset_b, offset_c;
    
    for(int i = 0; i < Order-1; ++i){
        for(int j = n; j < Order-1; ++j){
            BaseInt inter[Size];
            #pragma unroll
            for(std::size_t k=0 ; k < Size ;++k)
                inter[k] = 0;
            
            offset_a = (offset_block_ai+n-j)*Size;
            offset_b = (offset_block_bj-n)*Size;
            offset_c = ((offset_block_ai+offset_block_bj)*2-2*Order+2-j)*Size;
            
            single_multiplication_device<BaseInt,Size>(&a[offset_a],&b[offset_b],&inter[0]);
            kernels_addition_classic<BaseInt,Size>(&c[offset_c],&inter[0]);
        }    
        n++;
    }         
}

template <typename BaseInt, std::size_t Size>
void algo_block_algo(int i, int j, unsigned int Order, BaseInt const* a, BaseInt const* b, BaseInt *c)
{
    algo_triangle_up<BaseInt,Size>(i,j,Order,a,b,c);
    algo_diag<BaseInt,Size>(i,j,Order,a,b,c);
    algo_triangle_down<BaseInt,Size>(i,j,Order,a,b,c);
}


/**
* VLI_GPU_VECTOR functions
*/    
template  <typename BaseInt, std::size_t Size>
__global__ void inner_prod_vector(unsigned int max_order, std::size_t vector_size, BaseInt const* v1, BaseInt const* v2, BaseInt* res)
{
    unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x // get poly one by one
    const std::size_t size_multiplicant = Size*max_order*max_order;
    const std::size_t size_product = Size*2*max_order*2*max_order;
    if(xIndex < vector_size){
        //multiplication between polynomial
        std::size_t offset_m = xIndex*size_multiplicant;
        std::size_t offset_p = xIndex*size_product;
        polynome_polynome_multiplication_device<BaseInt,Size>(max_order,&v1[offset_m],&v2[offset_m],&res[offset_p]); 
    }
}
    
    
template  <typename BaseInt, std::size_t Size>
__global__ void reduction_polynome(unsigned int max_order, std::size_t vector_size, BaseInt* v1)
{ 
    std::size_t size_poly = Size*max_order*max_order;
    for(std::size_t i=1 ; i < vector_size ; ++i){
        for(std::size_t j=0 ; j < max_order*max_order; ++j){ //additional loop
            std::size_t offset0 = j*Size;
            std::size_t offset1 = i*size_poly+j*Size;
            kernels_addition_classic<BaseInt,Size>(&v1[offset0],&v1[offset1]);
        }
    }
}
    
/**
* New algo 
*/
    
template <typename BaseInt, std::size_t Size>
__global__ void inner_prod_vector_blocks(unsigned int Order, std::size_t vector_size, BaseInt const* A, BaseInt const* B, BaseInt* C)
{
    // remove the loops  
    unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x // get poly one by one
    unsigned int yindex = threadIdx.y; // thread for the triangle/diag decomposition
    const std::size_t size_multiplicant = Size*Order*Order;
    const std::size_t size_product = Size*2*Order*2*Order;
    //multiplication between polynomial
    std::size_t offset_m = xIndex*size_multiplicant;
    std::size_t offset_p = xIndex*size_product;
    if(xIndex < vector_size){
        // first pass, half top right corner, 
        for(int j=0; j<=yindex; ++j)
            algo_block_algo<BaseInt, Size>(j,yindex-j,Order,&A[offset_m],&B[offset_m],&C[offset_p]);
        
        //second pass, half bottom left corner
        for(int j=yindex+1; j<Order; ++j)
            algo_block_algo<BaseInt, Size>(j,Order-j+yindex,Order,&A[offset_m],&B[offset_m],&C[offset_p]);          
    }
}
            
     
/**
  * The C++ functions that call the kernels
  */
    
template <typename BaseInt, std::size_t Size>
void inner_product_vector(unsigned int Order, std::size_t vector_size, BaseInt const* A, BaseInt const* B, BaseInt* C, std::size_t threads_per_block) 
{
    std::size_t blocks_per_grid = vector_size/threads_per_block+1;
    dim3 dimgrid(blocks_per_grid,1,1);
    dim3 dimblock(threads_per_block,1,1);
    inner_prod_vector<BaseInt, Size> <<< dimgrid, dimblock >>>(Order, vector_size, A, B, C);
}
 
template <typename BaseInt, std::size_t Size>
void inner_product_vector_blocks(unsigned int Order, std::size_t vector_size, BaseInt const* A, BaseInt const* B, BaseInt *C)
{
   // std::size_t threads_per_block=64;
    
    dim3 dimgrid(vector_size,1,1);
    dim3 dimblock(1,Order,1);
    inner_prod_vector_blocks<BaseInt,Size><<<dimgrid,dimblock>>>(Order,vector_size,A,B,C);       
}
    
template <typename BaseInt, std::size_t Size>
void vector_reduction_inplace(unsigned int max_order, std::size_t vector_size, BaseInt* A)
{
    //the reduction should be // if executed on one smp
    dim3 dimgrid(1,1,1);
    dim3 dimblock(1,1,1);
    reduction_polynome<BaseInt, Size> <<< dimgrid, dimblock >>>(max_order, vector_size, A);
}        

#define VLI_IMPLEMENT_GPU_FUNCTIONS(TYPE, VLI_SIZE) \
    void inner_product_vector(vli_size_tag<VLI_SIZE>, unsigned int max_order, std::size_t vector_size, TYPE const* A, TYPE const* B, TYPE* C, std::size_t threads_per_block) \
        {inner_product_vector<TYPE,VLI_SIZE>(max_order,vector_size,A,B,C,threads_per_block);} \
    void vector_reduction_inplace(vli_size_tag<VLI_SIZE>, unsigned int max_order, std::size_t vector_size, TYPE* A) \
        {vector_reduction_inplace<TYPE,VLI_SIZE>(max_order,vector_size,A);} \
    void inner_product_vector_blocks(vli_size_tag<VLI_SIZE>, unsigned int Order, std::size_t vector_size, TYPE const* A, TYPE const* B, TYPE* C) \
        {inner_product_vector_blocks<TYPE,VLI_SIZE>(Order, vector_size, A, B, C);} 

    
#define VLI_IMPLEMENT_GPU_FUNCTIONS_FOR(r, data, BASEINT_SIZE_PAIR) \
    VLI_IMPLEMENT_GPU_FUNCTIONS( BOOST_PP_TUPLE_ELEM(2,0,BASEINT_SIZE_PAIR), BOOST_PP_TUPLE_ELEM(2,1,BASEINT_SIZE_PAIR) )

BOOST_PP_SEQ_FOR_EACH(VLI_IMPLEMENT_GPU_FUNCTIONS_FOR, _, VLI_COMPILE_BASEINT_SIZE_PAIRS_SEQ)

#undef VLI_IMPLEMENT_GPU_FUNCTIONS_FOR
#undef VLI_IMPLEMENT_GPU_FUNCTIONS

} // namespace detail
} // namespace vli

