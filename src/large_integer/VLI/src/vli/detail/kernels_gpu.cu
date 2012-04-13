//
//  kernels_gpu_interface.h
//  vli
//
//  Created by Timothée Ewart on 22/08/11.
//  Copyright 2011 Université de Genève. All rights reserved.
//

#include "vli/detail/kernels_gpu_asm.hpp" 
#include "vli/detail/kernels_gpu.h"

#include <cstdio>
//#include "utils/cuPrintf.cu"

namespace vli {
    namespace detail {
    
// functions for the blocks algos
/*
    template <typename BaseInt, std::size_t Size>
    __device__ void algo_triangle_up(int block_ai, int block_bj,unsigned int Order, BaseInt const* a,  BaseInt const* b, BaseInt* c);

    template <typename BaseInt, std::size_t Size>
    __device__ void algo_triangle_down(int block_ai, int block_bj,unsigned int Order, BaseInt const* a,  BaseInt const* b, BaseInt* c);
        
    template <typename BaseInt, std::size_t Size>
    __device__ void algo_diag(int block_ai, int block_bj,unsigned int Order, BaseInt const* a, BaseInt const* b, BaseInt* c);
    
    template <typename BaseInt, std::size_t Size>
    __device__ void algo_block_algo(int i, int j,unsigned int Order, BaseInt const* a, BaseInt const* b, BaseInt* c);
*/
// functions for the diags algo
/*
    template <typename BaseInt, std::size_t Size>
    __device__ void algo_diag_shared(unsigned int i,unsigned int Order, BaseInt const* a,  BaseInt const* b, BaseInt* c);
*/
/**
* a kind of hook function with a little bit of arithmetic in case of signed int (multiplication)
*/

/**
* VLI_GPU functions
*/

/** n threads block algo, note still with truncated multiplication **/        
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
             
       //     single_multiplication_device<BaseInt,Size>(&a[offset_a],&b[offset_b],inter);
       //     kernels_addition_classic<BaseInt,Size>(&c[offset_c],inter);
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

     //   single_multiplication_device<BaseInt,Size>(&a[offset_a],&b[offset_b],inter);
     //   kernels_addition_classic<BaseInt,Size>(&c[offset_c],inter);
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
            
        //    single_multiplication_device<BaseInt,Size>(&a[offset_a],&b[offset_b],&inter[0]);
        //    kernels_addition_classic<BaseInt,Size>(&c[offset_c],&inter[0]);
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

/** n*n threads diag algo **/
/*
template <typename BaseInt, std::size_t Size>
void algo_diag_shared(unsigned int threadid, unsigned int Order, BaseInt const* a, BaseInt const* b, BaseInt *c)
{
    int qa,ra,qb,rb,pos; // find all indexes
    int offset_a, offset_b, offset_c;
    int j = Order*Order-1;
    unsigned int oldthreadid(threadid);

    __shared__ BaseInt sc[2*Size*121];
    __shared__ BaseInt s_inter[2*Size*121];
    __shared__ BaseInt sa[Size*121]; // 121 = ORDER*ORDER
    __shared__ BaseInt sb[Size*121]; // 121 = ORDER*ORDER
 
    for(int i(0); i <= threadid; i++){

        qa = i/Order;
        ra = i%Order;
        qb = (threadid-i)/Order;
        rb = (threadid-i)%Order;
        pos = 2*(qa+qb)*Order + (ra+rb);

        offset_a = (threadid-i)*Size;
        offset_b = i*Size;
        offset_c = pos*Size*2; // 2 because non truncated

        // we load in a faster memory
        #pragma unroll 
        for(std::size_t k=0 ; k < Size ;++k){
            sa[threadid*Size+k] = a[offset_a+k];
            sb[threadid*Size+k] = b[offset_b+k];
        }

        #pragma unroll 
        for(std::size_t k=0 ; k < 2*Size ;++k){ // 2 because non truncated
            s_inter[2*threadid*Size+k] = 0;
            sc[2*threadid*Size+k] = c[offset_c+k]; 
        }

 //       single_multiplication_device<BaseInt,Size>(&sa[threadid*Size],&sb[threadid*Size],&s_inter[2*threadid*Size]);        
 //       kernels_addition_classic<BaseInt,2*Size>(&sc[2*threadid*Size],&s_inter[2*threadid*Size]); // 2 because non truncated

        #pragma unroll 
        for(std::size_t k=0 ; k < 2*Size ;++k) // 2 because non truncated
           c[offset_c+k] = sc[2*threadid*Size+k]; 
      
        __syncthreads();
    }
    
    threadid  = Order*Order - threadid ; // back flip of the thread

    for(int i(Order*Order-threadid+1); i < Order*Order; i++){
        #pragma unroll
        for(std::size_t k=0 ; k < 2*Size ;++k) // 2 because non truncated
            s_inter[k] = 0;

        qa = i/Order;
        ra = i%Order;
        qb = j/Order;
        rb = j%Order;
        pos = 2*(qa+qb)*Order + (ra+rb);

        offset_a = j*Size;
        offset_b = i*Size;
        offset_c = pos*Size*2; // 2 because non truncated
    
        // we load in a faster memory
        #pragma unroll 
        for(std::size_t k=0 ; k < Size ;++k){
            sa[oldthreadid*Size+k] = a[offset_a+k];
            sb[oldthreadid*Size+k] = b[offset_b+k];
        }

        #pragma unroll 
        for(std::size_t k=0 ; k < 2*Size ;++k){ // 2 because non truncated
            s_inter[2*oldthreadid*Size+k] = 0;
            sc[2*oldthreadid*Size+k] = c[offset_c+k]; 
        }

 //       single_multiplication_device<BaseInt,Size>(&sa[oldthreadid*Size],&sb[oldthreadid*Size],&s_inter[2*oldthreadid*Size]);        
  //      kernels_addition_classic<BaseInt,2*Size>(&sc[2*oldthreadid*Size],&s_inter[2*oldthreadid*Size]); // 2 because non truncated

        #pragma unroll 
        for(std::size_t k=0 ; k < 2*Size ;++k) // 2 because non truncated
           c[offset_c+k] = sc[2*oldthreadid*Size+k]; 

        __syncthreads();
        j--;        
    }    

}
*/
/**
* VLI_GPU_VECTOR functions
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
* New algo based on diagonal decomposition 
*/

template <typename BaseInt, std::size_t Size>
__global__ void inner_prod_vector_diag(unsigned int Order, std::size_t vector_size, BaseInt const* A, BaseInt const* B, BaseInt* C)
{
    // remove the loops  
    unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x // get poly one by one
    unsigned int yIndex = threadIdx.y; // thread for the triangle/diag decomposition
    const std::size_t size_multiplicant = Size*Order*Order;
    const std::size_t size_product = 2*Size*2*Order*2*Order;
    //multiplication between polynomial
    std::size_t offset_m = xIndex*size_multiplicant;
    std::size_t offset_p = xIndex*size_product;
    if(xIndex < vector_size){
      // algo_diag_shared<BaseInt,Size>(yIndex,Order,&A[offset_m],&B[offset_m],&C[offset_p]);
      //  algo_diag<BaseInt,Size>(yIndex,Order,&A[offset_m],&B[offset_m],&C[offset_p]);
    }
}

// C - on cpu we use unsigned long int 
template  <typename BaseInt, std::size_t Size, std::size_t ThreadsPerBlock>
__global__ void reduction_polynome(unsigned int Order, std::size_t VectorSize, BaseInt* v)
{ 
// C - tody this kernel 0.02 s ....
    unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x // get poly one by one
    unsigned int xInter = 2*Size*threadIdx.x;  
    unsigned int offset_poly  = 2*Size*2*Order*2*Order;
    unsigned int offset_coeff = 2*Size*xIndex; 

    __shared__ unsigned long int sv[2*Size*ThreadsPerBlock];

    if(xIndex < 2*Order*2*Order){
        for(int j(0); j < 2*Size ; ++j) 
            sv[xInter+j] = v[offset_coeff+j];
/*
        for(int i(1); i < VectorSize; ++i)
            add384_384_gpu(&sv[xInter], &v[offset_coeff+i*offset_poly]);
*/
        for(int j(0); j < 2*Size ; ++j) 
            v[offset_coeff+j] = sv[xInter+j];
    }
}

template  <typename BaseInt, std::size_t Size>
__global__ void addition(BaseInt* x, BaseInt const* y){ 
    add384_384_gpu(x,y);
}


template  <typename BaseInt, std::size_t Size>
__global__ void multiplication(BaseInt* x, BaseInt const* y, BaseInt const* z){ 
    mul384_384_gpu(x,y,z);
}


/**
  * The C++ functions that call the kernels
*/
template <typename BaseInt, std::size_t Size>
void inner_product_vector(unsigned int Order, std::size_t VectorSize, BaseInt const* A, BaseInt const* B, BaseInt *C) {

  std::size_t ThreadsPerBlock = 32;

  dim3 dimgrid(VectorSize,1,1);
  dim3 dimblock(1,Order*Order,1);
//  inner_prod_vector_blocks<BaseInt,Size><<<dimgrid,dimblock>>>(Order,VectorSize,A,B,C);      // nthreads version truncated multiplication 
//  inner_prod_vector_diag<BaseInt,Size><<<dimgrid,dimblock>>>(Order,VectorSize,A,B,C);          // nthreads*nthreads version non truncated 
  //change the grid size
  dimgrid.x  = 2*Order*2*Order/ThreadsPerBlock+1;
  dimblock.x = ThreadsPerBlock;
  dimblock.y = 1;
  reduction_polynome<BaseInt,Size, 32> <<<dimgrid, dimblock>>>(Order, VectorSize,C);
}

// just for testing my asm kernel
// test addition
template <typename BaseInt, std::size_t Size>
void addition_gpu(BaseInt* A, BaseInt const* B) {
  dim3 dimgrid(1,1,1);
  dim3 dimblock(1,1,1);
  addition<BaseInt,Size> <<<dimgrid, dimblock>>>(A, B);
}
//test multiplication
template <typename BaseInt, std::size_t Size>
void multiplication_gpu(BaseInt* A, BaseInt const* B, BaseInt const* C) {
  dim3 dimgrid(1,1,1);
  dim3 dimblock(1,1,1);
  multiplication<BaseInt,Size> <<<dimgrid, dimblock>>>(A, B, C);
}
// the worst cast of my life
#define VLI_IMPLEMENT_GPU_FUNCTIONS(TYPE, VLI_SIZE) \
    void inner_product_vector(vli_size_tag<VLI_SIZE>, unsigned int max_order, std::size_t vector_size, TYPE const* A, TYPE const* B, TYPE* C) \
        {inner_product_vector<TYPE,VLI_SIZE>(max_order,vector_size,A,B,C);} \
    void addition_gpu(vli_size_tag<2*VLI_SIZE>, TYPE* A, TYPE const* B) \
        {addition_gpu<unsigned int,2*VLI_SIZE>((unsigned int*)A,(unsigned int*)(void*)B);} \
    void multiplication_gpu(vli_size_tag<VLI_SIZE>, TYPE* A, TYPE const* B, TYPE const* C) \
        {multiplication_gpu<unsigned int,2*VLI_SIZE>((unsigned int*)A,(unsigned int*)(void*)B,(unsigned int*)(void*)C);} \

#define VLI_IMPLEMENT_GPU_FUNCTIONS_FOR(r, data, BASEINT_SIZE_PAIR) \
    VLI_IMPLEMENT_GPU_FUNCTIONS( BOOST_PP_TUPLE_ELEM(2,0,BASEINT_SIZE_PAIR), BOOST_PP_TUPLE_ELEM(2,1,BASEINT_SIZE_PAIR) )

BOOST_PP_SEQ_FOR_EACH(VLI_IMPLEMENT_GPU_FUNCTIONS_FOR, _, VLI_COMPILE_BASEINT_SIZE_PAIRS_SEQ)

#undef VLI_IMPLEMENT_GPU_FUNCTIONS_FOR
#undef VLI_IMPLEMENT_GPU_FUNCTIONS

} // namespace detail
} // namespace vli

