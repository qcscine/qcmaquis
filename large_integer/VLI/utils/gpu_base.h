//
//  gpu_p.h
//  XCODE_VLI
//
//  Created by Tim Ewart on 07.07.11.
//  Copyright 2011 University of Geneva. All rights reserved.
//


#ifndef __VLI_GPU_POINTER_ARRAY_VECTOR__
#define __VLI_GPU_POINTER_ARRAY_VECTOR__

#include <iostream>
#include <boost/static_assert.hpp>
#include "gpu/GpuManager.h"

namespace vli
{   
    template<class T>
    class gpu_pointer
    {
    public:
        
        gpu_pointer():data_(NULL){};
        
        ~gpu_pointer(){}
        
        
        T* p(){
            return data_;
        }
        
        T* const p() const{
            return data_;
        }
        
    protected: // for the swap, to change ?
        T* data_;
    };

    template<class BaseInt, std::size_t Size>
    class gpu_array : public gpu_pointer<BaseInt>
    {
    public:
        gpu_array(){ // here, Size is template
		    gpu::cu_check_error(cudaMalloc((void**)&(this->data_), Size*sizeof(BaseInt)), __LINE__);			
	    	gpu::cu_check_error(cudaMemset((void*)(this->data_), 0, Size*sizeof(BaseInt)), __LINE__);			
        }
                
        ~gpu_array(){
             gpu::cu_check_error(cudaFree(this->data_), __LINE__);   
        }
    };

    template<class BaseInt>
    struct gpu_vector : public gpu_pointer<BaseInt>
    {
    public:
        gpu_vector(std::size_t Size){ // here, Size is an argument
		    gpu::cu_check_error(cudaMalloc((void**)&(this->data_), Size*sizeof(BaseInt)), __LINE__);			
	    	gpu::cu_check_error(cudaMemset((void*)(this->data_), 0, Size*sizeof(BaseInt)), __LINE__);			            
        }
        
        ~gpu_vector(){
            gpu::cu_check_error(cudaFree(this->p()), __LINE__);   
        }
        
        void resize(std::size_t newsize){
            BaseInt* temp;
            gpu::cu_check_error(cudaMalloc((void**)(&temp),newsize*sizeof(BaseInt)), __LINE__);
            gpu::cu_check_error(cudaMemset((void*)(temp),0, newsize*sizeof(BaseInt)), __LINE__);	
            gpu::cu_check_error(cudaMemcpy((void*)(temp),this->data_, newsize*sizeof(BaseInt),cudaMemcpyDeviceToDevice), __LINE__); 
            std::swap(temp,this->data_);
            gpu::cu_check_error(cudaFree(temp), __LINE__);
        }
        
        const std::size_t size() const {
            return Size_;
        }
        
    private:
        std::size_t Size_; // number of elements
    };
    
}

#endif