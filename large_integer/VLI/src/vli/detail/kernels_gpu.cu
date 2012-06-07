//
//  kernels_gpu_interface.h
//  vli
//
//  Created by Timothée Ewart on 22/08/11.
//  Copyright 2011 Université de Genève. All rights reserved.
//

#include "vli/detail/numeric.h"
#include "vli/detail/kernels_gpu.h"
#include "vli/detail/kernels_gpu_asm.hpp" 
#include "vli/detail/diag_algo_ip_gpu.h"
#include <stdio.h>
namespace vli {
    namespace detail {
  
    /** 
    * Algo based on diagonal decomposition 
    */
    template <typename BaseInt, std::size_t Size, unsigned int Order>
    __global__ void inner_prod_vector_diag(std::size_t vector_size, BaseInt const* A, BaseInt const* B, BaseInt* C) {
        unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x // get poly one by one
        unsigned int yIndex = threadIdx.y; // thread for the triangle/diag decomposition
        const std::size_t size_multiplicant = Size*Order*Order;
        const std::size_t size_product = 2*Size*2*Order*2*Order;
        //multiplication between polynomial
        std::size_t offset_m = xIndex*size_multiplicant;
        std::size_t offset_p = xIndex*size_product;
    
        if(xIndex < vector_size)
            if(yIndex < Order*Order)
                diag_algo<BaseInt,Size, Order>(yIndex,&A[offset_m],&B[offset_m],&C[offset_p]);

    } 

    template  <typename BaseInt, std::size_t Size, std::size_t ThreadsPerBlock>
     __global__ void reduction_polynome(unsigned int Order, std::size_t VectorSize, BaseInt* v) { 
        unsigned int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x // get poly one by one
        unsigned int xInter = 2*Size*threadIdx.x;  
        unsigned int offset_poly  = 2*Size*2*Order*2*Order;
        unsigned int offset_coeff = 2*Size*xIndex; 
   
        __shared__ unsigned int sv[2*Size*ThreadsPerBlock];
   
        if(xIndex < 2*Order*2*Order){
            for(int j(0); j < 2*Size ; ++j) 
                sv[xInter+j] = v[offset_coeff+j];
    
            for(int i(1); i < VectorSize; ++i)
                add384_384_gpu(&sv[xInter], &v[offset_coeff+i*offset_poly]);
   
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
   
    template <typename BaseInt, std::size_t Size, unsigned int Order>
    void inner_product_vector(std::size_t VectorSize, BaseInt const* A, BaseInt const* B, BaseInt *C) {
   
      std::size_t ThreadsPerBlock = 32;
      dim3 dimgrid(VectorSize,1,1);
      dim3 dimblock(1,Order*Order,1);
      inner_prod_vector_diag<BaseInt,Size,Order><<<dimgrid,dimblock>>>(VectorSize,A,B,C);       
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

    #define VLI_IMPLEMENT_GPU_FUNCTIONS(TYPE, VLI_SIZE, POLY_ORDER) \
        void inner_product_vector(vli_size_tag<VLI_SIZE, POLY_ORDER>, std::size_t vector_size, TYPE const* A, TYPE const* B, TYPE* C) \
            {inner_product_vector<unsigned int, 2*VLI_SIZE, POLY_ORDER>(vector_size, const_cast<unsigned int*>(reinterpret_cast<unsigned int const*>(A)), const_cast<unsigned int*>(reinterpret_cast<unsigned int const*>(B)),(unsigned int*)C);} \
        void addition_gpu(vli_size_tag<2*VLI_SIZE,POLY_ORDER>, TYPE* A, TYPE const* B) \
            {addition_gpu<unsigned int,2*VLI_SIZE>((unsigned int*)A,const_cast<unsigned int*>(reinterpret_cast<unsigned int const*>(B)));} \
        void multiplication_gpu(vli_size_tag<VLI_SIZE,POLY_ORDER>, TYPE* A, TYPE const* B, TYPE const* C) \
            {multiplication_gpu<unsigned int,2*VLI_SIZE>((unsigned int*)A, const_cast<unsigned int*>(reinterpret_cast<unsigned int const*>(B)),(const_cast<unsigned int*>(reinterpret_cast<unsigned int const*>(C))));} 
   
    #define VLI_IMPLEMENT_GPU_FUNCTIONS_FOR(r, data, BASEINT_SIZE_ORDER_TUPLE) \
        VLI_IMPLEMENT_GPU_FUNCTIONS( BOOST_PP_TUPLE_ELEM(3,0,BASEINT_SIZE_ORDER_TUPLE), BOOST_PP_TUPLE_ELEM(3,1,BASEINT_SIZE_ORDER_TUPLE), BOOST_PP_TUPLE_ELEM(3,2,BASEINT_SIZE_ORDER_TUPLE) )
   
    BOOST_PP_SEQ_FOR_EACH(VLI_IMPLEMENT_GPU_FUNCTIONS_FOR, _, VLI_COMPILE_BASEINT_SIZE_ORDER_TUPLE_SEQ)
   
    #undef VLI_IMPLEMENT_GPU_FUNCTIONS_FOR
    #undef VLI_IMPLEMENT_GPU_FUNCTIONS

} // namespace detail
} // namespace vli

