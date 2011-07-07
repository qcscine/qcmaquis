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
        
        gpu_array(gpu_array const& a)
        {
		    gpu::cu_check_error(cudaMalloc((void**)&(this->data_), Size*sizeof(BaseInt)), __LINE__);
            gpu::cu_check_error(cudaMemcpy((void*)this->data_,(void*)a.data_,Size*sizeof(BaseInt), cudaMemcpyDeviceToDevice), __LINE__);
        } 
        
        ~gpu_array(){
             gpu::cu_check_error(cudaFree(this->data_), __LINE__);   
        }

        gpu_array& operator = (gpu_array a)
        {
            swap(*this,a);
            return *this;
        }

        friend void swap(gpu_array& a1, gpu_array& a2)
        {
            std::swap(a1.data_,a2.data_);
        }
    };

    template<class BaseInt>
    struct gpu_vector : public gpu_pointer<BaseInt>
    {
    public:
        gpu_vector(std::size_t size = 8)
        :size_(size){ // here, Size is an argument
            gpu::cu_check_error(cudaMalloc((void**)&(this->data_), size*sizeof(BaseInt)), __LINE__);
            gpu::cu_check_error(cudaMemset((void*)this->data_,0, size*sizeof(BaseInt)), __LINE__);   
        }
        
        gpu_vector(gpu_vector const& a)
        :size_(a.size_)
        {
            gpu::cu_check_error(cudaMalloc((void**)&(this->data_), size_*sizeof(BaseInt)), __LINE__);
            gpu::cu_check_error(cudaMemcpy((void*)this->data_,(void*)a.data_,a.size_*sizeof(BaseInt), cudaMemcpyDeviceToDevice), __LINE__);
        }
        
        ~gpu_vector(){
            gpu::cu_check_error(cudaFree(this->data_), __LINE__);
        }
        
        void resize(std::size_t newsize){
            BaseInt* temp;
            gpu::cu_check_error(cudaMalloc((void**)(&temp), newsize*sizeof(BaseInt)), __LINE__);
            if(newsize > size_)
            {
                gpu::cu_check_error(cudaMemcpy((void*)(temp), this->data_, size_*sizeof(BaseInt),cudaMemcpyDeviceToDevice), __LINE__);
                gpu::cu_check_error(cudaMemset((void*)(temp+size_),0, (newsize-size_)*sizeof(BaseInt)), __LINE__);
            }
            else
            {
                gpu::cu_check_error(cudaMemcpy((void*)(temp),this->data_, newsize*sizeof(BaseInt),cudaMemcpyDeviceToDevice), __LINE__);
            }
            std::swap(temp,this->data_);
            gpu::cu_check_error(cudaFree(temp), __LINE__);
        }
        
        gpu_vector& operator = (gpu_vector a)
        {
            swap(*this,a);
            return *this;
        }
        
        friend void swap(gpu_vector& a1, gpu_vector& a2)
        {
            using std::swap;
            swap(a1.size_,a2.size_);
            swap(a1.data_,a2.data_);
        }
        
        std::size_t size() const {
            return size_;
        }
 
    private:
        std::size_t size_; // number of elements
    };
    
}

#endif
