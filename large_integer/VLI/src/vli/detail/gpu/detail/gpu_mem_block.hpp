/*
*Very Large Integer Library, License - Version 1.0 - May 3rd, 2012
*
*Timothee Ewart - University of Geneva, 
*Andreas Hehn - Swiss Federal Institute of technology Zurich.
*Maxim Milakov - NVIDIA
*
*Permission is hereby granted, free of charge, to any person or organization
*obtaining a copy of the software and accompanying documentation covered by
*this license (the "Software") to use, reproduce, display, distribute,
*execute, and transmit the Software, and to prepare derivative works of the
*Software, and to permit third-parties to whom the Software is furnished to
*do so, all subject to the following:
*
*The copyright notices in the Software and this entire statement, including
*the above license grant, this restriction and the following disclaimer,
*must be included in all copies of the Software, in whole or in part, and
*all derivative works of the Software, unless such copies or derivative
*works are solely in the form of machine-executable object code generated by
*a source language processor.
*
*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
*SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
*FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
*ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
*DEALINGS IN THE SOFTWARE.
*/
// to remove when size ok
#ifndef GPU_MEM_BLOCK_HPP
#define GPU_MEM_BLOCK_HPP

#include <boost/mpl/assert.hpp>
#include "vli/detail/gpu/detail/gpu_mem_block.h" // memory

namespace vli{
namespace detail {

    gpu_memblock::gpu_memblock()
    : block_size_(0), V1Data_(0), V2Data_(0), VinterData_(0), PoutData_(0) {
    }

    std::size_t const& gpu_memblock::BlockSize() const {
        return block_size_;
    }; 

    template <std::size_t NumBits, class MaxOrder, int NumVars>
    struct resize_helper{
    };
    // A free function, syntax light compared nested into the class
    // max order each specialization 
    template <std::size_t NumBits, int Order, int NumVars>
    struct resize_helper<NumBits, max_order_each<Order>, NumVars>{
        static void resize(gpu_memblock const& pgm,  std::size_t vector_size){

        std::size_t req_size = vector_size * num_words<NumBits>::value * stride<0,NumVars,Order>::value * stride<1,NumVars,Order>::value * stride<2,NumVars,Order>::value * stride<3,NumVars,Order>::value;

        if( req_size > pgm.BlockSize() ) {
            pgm.block_size_  = req_size;
            if(pgm.V1Data_ != 0 )
                gpu::cu_check_error(cudaFree((void*)pgm.V1Data_),__LINE__);
            if(pgm.V2Data_ != 0 )
                gpu::cu_check_error(cudaFree((void*)pgm.V2Data_),__LINE__);
            if(pgm.VinterData_ != 0)
                gpu::cu_check_error(cudaFree((void*)pgm.VinterData_),__LINE__);
            if(pgm.PoutData_ != 0)
                gpu::cu_check_error(cudaFree((void*)pgm.PoutData_),__LINE__);
            gpu::cu_check_error(cudaMalloc((void**)&(pgm.V1Data_), req_size*sizeof(boost::uint32_t)),__LINE__); //input 1
            gpu::cu_check_error(cudaMalloc((void**)&(pgm.V2Data_), req_size*sizeof(boost::uint32_t)),__LINE__); //input 2
            gpu::cu_check_error(cudaMalloc((void**)&(pgm.VinterData_), vector_size *  2*num_words<NumBits>::value 
                                                     * result_stride<0, NumVars, Order>::value * result_stride<1, NumVars, Order>::value * result_stride<2, NumVars, Order>::value * result_stride<3, NumVars, Order>::value * sizeof(boost::uint32_t)),__LINE__); 
            gpu::cu_check_error(cudaMalloc((void**)&(pgm.PoutData_),                  2*num_words<NumBits>::value 
                                                     * result_stride<0, NumVars, Order>::value * result_stride<1, NumVars, Order>::value * result_stride<2, NumVars, Order>::value * result_stride<3, NumVars, Order>::value * sizeof(boost::uint32_t)),__LINE__);
        } // end if
        gpu::cu_check_error(cudaMemset((void*)(pgm.PoutData_),0,2*num_words<NumBits>::value
                                               * result_stride<0, NumVars, Order>::value * result_stride<1, NumVars, Order>::value * result_stride<2, NumVars, Order>::value * result_stride<3, NumVars, Order>::value * sizeof(boost::uint32_t)),__LINE__);       

        } // end function
    }; // end struct

    // max order combined specialization 
    template <std::size_t NumBits, int Order, int NumVars>
    struct resize_helper<NumBits, max_order_combined<Order>, NumVars>{
        static void resize(gpu_memblock const& pgm,  std::size_t vector_size){

        std::size_t req_size = vector_size * num_words<NumBits>::value * vli::detail::max_order_combined_helpers::size<NumVars+1, Order>::value;

        if( req_size > pgm.BlockSize() ) {
            pgm.block_size_  = req_size;
            if(pgm.V1Data_ != 0 )
                gpu::cu_check_error(cudaFree((void*)pgm.V1Data_),__LINE__);
            if(pgm.V2Data_ != 0 )
                gpu::cu_check_error(cudaFree((void*)pgm.V2Data_),__LINE__);
            if(pgm.VinterData_ != 0)
                gpu::cu_check_error(cudaFree((void*)pgm.VinterData_),__LINE__);
            if(pgm.PoutData_ != 0)
                gpu::cu_check_error(cudaFree((void*)pgm.PoutData_),__LINE__);

            gpu::cu_check_error(cudaMalloc((void**)&(pgm.V1Data_), req_size*sizeof(boost::uint32_t)),__LINE__); //input 1
            gpu::cu_check_error(cudaMalloc((void**)&(pgm.V2Data_), req_size*sizeof(boost::uint32_t)),__LINE__); //input 2
            gpu::cu_check_error(cudaMalloc((void**)&(pgm.VinterData_), vector_size *  2*num_words<NumBits>::value * vli::detail::max_order_combined_helpers::size<NumVars+1, 2*Order>::value*sizeof(boost::uint32_t)),__LINE__); 
            gpu::cu_check_error(cudaMalloc((void**)&(pgm.PoutData_),                  2*num_words<NumBits>::value * vli::detail::max_order_combined_helpers::size<NumVars+1, 2*Order>::value*sizeof(boost::uint32_t)),__LINE__);

            } // end if
        gpu::cu_check_error(cudaMemset((void*)(pgm.PoutData_),0,2*num_words<NumBits>::value * vli::detail::max_order_combined_helpers::size<NumVars+1, 2*Order>::value*sizeof(boost::uint32_t)),__LINE__);
        } // end fonction
    }; //end struct


    template <std::size_t NumBits, class MaxOrder, int NumVars>
    struct memory_transfer_helper;
    
    // max order each specialization 
    template <std::size_t NumBits, int Order, int NumVars>
    struct memory_transfer_helper<NumBits, max_order_each<Order>, NumVars>{
         static void transfer_up(gpu_memblock const& pgm, boost::uint32_t const* pData1, boost::uint32_t const* pData2,  std::size_t VectorSize){
            BOOST_STATIC_ASSERT( NumVars <= 2);
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V1Data_,(void*)pData1,VectorSize*stride<0,NumVars,Order>::value*stride<1,NumVars,Order>::value*stride<2,NumVars,Order>::value*stride<3,NumVars,Order>::value
                                *num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V2Data_,(void*)pData2,VectorSize*stride<0,NumVars,Order>::value*stride<1,NumVars,Order>::value*stride<2,NumVars,Order>::value*stride<3,NumVars,Order>::value
                                *num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);
         }
    };

    template <std::size_t NumBits, int Order>
    struct memory_transfer_helper<NumBits, max_order_each<Order>, 3 >{
         static void transfer_up(gpu_memblock const& pgm, boost::uint32_t const* pData1, boost::uint32_t const* pData2,  std::size_t VectorSize){
            
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V1Data_,(void*)pData1,VectorSize*stride<0,3,Order>::value*stride<1,3,Order>::value*stride<2,3,Order>::value*stride<3,3,Order>::value
                                *num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V2Data_,(void*)pData2,VectorSize*stride<0,3,Order>::value*stride<1,3,Order>::value*stride<2,3,Order>::value*stride<3,3,Order>::value
                                *num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);

            //only the second poly is cashed into the texture memory
            gpu::cu_check_error(cudaBindTexture(0,tex_reference_2,(void*)pgm.V2Data_,VectorSize*stride<0,3,Order>::value*stride<1,3,Order>::value*stride<2,3,Order>::value*stride<3,3,Order>::value
                                *num_words<NumBits>::value*sizeof(boost::uint32_t)),__LINE__);

         }
    };

    template <std::size_t NumBits, int Order>
    struct memory_transfer_helper<NumBits, max_order_each<Order>, 4 >{
         static void transfer_up(gpu_memblock const& pgm, boost::uint32_t const* pData1, boost::uint32_t const* pData2,  std::size_t VectorSize){
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V1Data_,(void*)pData1,VectorSize*stride<0,4,Order>::value*stride<1,4,Order>::value*stride<2,4,Order>::value*stride<3,4,Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V2Data_,(void*)pData2,VectorSize*stride<0,4,Order>::value*stride<1,4,Order>::value*stride<2,4,Order>::value*stride<3,4,Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);

            //two polys are cached, too large
            gpu::cu_check_error(cudaBindTexture(0,tex_reference_1,(void*)pgm.V1Data_,VectorSize*stride<0,4,Order>::value*stride<1,4,Order>::value*stride<2,4,Order>::value*stride<3,4,Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t)),__LINE__);
            gpu::cu_check_error(cudaBindTexture(0,tex_reference_2,(void*)pgm.V2Data_,VectorSize*stride<0,4,Order>::value*stride<1,4,Order>::value*stride<2,4,Order>::value*stride<3,4,Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t)),__LINE__);
         }
    };

    // max order combined  specialization 
    template <std::size_t NumBits, int Order, int NumVars>
    struct memory_transfer_helper<NumBits, max_order_combined<Order>, NumVars>{
         static void transfer_up(gpu_memblock const& pgm, boost::uint32_t const* pData1, boost::uint32_t const* pData2,  std::size_t VectorSize){
            BOOST_STATIC_ASSERT( NumVars <= 3);
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V1Data_,(void*)pData1,VectorSize*max_order_combined_helpers::size<NumVars+1, Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V2Data_,(void*)pData2,VectorSize*max_order_combined_helpers::size<NumVars+1, Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);
         }
    };

    template <std::size_t NumBits, int Order>
    struct memory_transfer_helper<NumBits, max_order_combined<Order>, 4>{
         static void transfer_up(gpu_memblock const& pgm, boost::uint32_t const* pData1, boost::uint32_t const* pData2,  std::size_t VectorSize){
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V1Data_,(void*)pData1,VectorSize*max_order_combined_helpers::size<4+1, Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);
  	    gpu::cu_check_error(cudaMemcpyAsync((void*)pgm.V2Data_,(void*)pData2,VectorSize*max_order_combined_helpers::size<4+1, Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t),cudaMemcpyHostToDevice),__LINE__);
            //only the second poly is cashed into the texture memory
            gpu::cu_check_error(cudaBindTexture(0,tex_reference_2,(void*)pgm.V2Data_,VectorSize*max_order_combined_helpers::size<4+1, Order>::value*num_words<NumBits>::value*sizeof(boost::uint32_t)),__LINE__);
         }
    };
 
    void UnbindTexture(){
        gpu::cu_check_error(cudaUnbindTexture(tex_reference_1),__LINE__);
        gpu::cu_check_error(cudaUnbindTexture(tex_reference_2),__LINE__);
    }

    } // end namespace detail
}// end namespace vli

#endif
