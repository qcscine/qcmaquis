/*
 *  vli.h
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef __VLI_GPU__
#define __VLI_GPU__


#include <iostream>
#include <boost/static_assert.hpp>
#include "function_hooks/vli_number_gpu_function_hooks.hpp"
#include "gpu/GpuManager.h"

namespace vli
{    
    template<class BaseInt, int Size>
    class vli_cpu;
    
	template<class BaseInt, int Size>
	class vli_gpu 
	{
	public:
	    typedef BaseInt         value_type;     // Data type to store parts of the very long integer (usually int)
        typedef std::size_t     size_type;      // Size type of the very long integers (number of parts)
        
        enum {size = Size};                     // Number of parts of the very long integer (eg. how many ints)
        
        /**
         proxy objects to access elements of the VLI
         */ 
        class proxy
        {
        public:
            proxy(BaseInt* p, size_type i)
            :data_(p), pos(i)
            {
            }
            
            operator vli_gpu() const
            {
                vli_gpu vli;
                gpu::cu_check_error( cudaMemcpy( vli.p(), data_, Size*sizeof(BaseInt) , cudaMemcpyDeviceToDevice), __LINE__);
                return vli;
            }
            
            proxy& operator= (BaseInt i)
            {
                gpu::cu_check_error(cudaMemcpy( (void*)(data_+pos),(void*)&i, sizeof(BaseInt), cudaMemcpyHostToDevice ), __LINE__); 					
                return *this;
            }
            
        private:
            size_type pos;
            BaseInt* data_;    
        };


        
		/**
		constructors 
		*/
        vli_gpu()
        {
		    gpu::cu_check_error(cudaMalloc((void**)&data_, Size*sizeof(BaseInt)), __LINE__);			
	    	gpu::cu_check_error(cudaMemset((void*)data_, 0, Size*sizeof(BaseInt)), __LINE__);			
        }

		explicit vli_gpu(int num)
        {
            gpu::cu_check_error(cudaMalloc((void**)&data_, Size*sizeof(BaseInt)), __LINE__);			
			gpu::cu_check_error(cudaMemset((void*)data_, 0, Size*sizeof(BaseInt)), __LINE__); //Everything is set to 0	
			gpu::cu_check_error(cudaMemset((void*)data_, num, sizeof(BaseInt)), __LINE__); //The first number is set to the num value			
        }

		vli_gpu(vli_gpu const& vli)
        {
			gpu::cu_check_error(cudaMalloc((void**)&data_, Size*sizeof(BaseInt)), __LINE__);
            gpu::cu_check_error(cudaMemcpy((void*)data_, vli.data_, Size*sizeof(BaseInt) , cudaMemcpyDeviceToDevice), __LINE__);
        }

		vli_gpu(vli_cpu<BaseInt, Size> const& vli)
        {
			gpu::cu_check_error(cudaMalloc((void**)&data_, Size*sizeof(BaseInt)  ), __LINE__);
			gpu::cu_check_error(cudaMemcpy((void*)data_, (void*)&vli[0], Size*sizeof(BaseInt), cudaMemcpyHostToDevice ), __LINE__); 		
        }
		
		/**
		destructors 
		*/
		~vli_gpu()
        {
            gpu::cu_check_error(cudaFree(data_), __LINE__);
        }

		friend void swap(vli_gpu& a, vli_gpu& b)
        {
            using std::swap;
            swap(a.data_, b.data_);
        }

		vli_gpu& operator = (vli_gpu vli)
        {
            // TODO perhaps the copy-swap implementation is not optimal on the GPU, to check
            swap(*this,vli);
            return *this;
        }

		void copy_vli_to_cpu(vli::vli_cpu<BaseInt, Size>& vli) const
        {
 			gpu::cu_check_error(cudaMemcpy( (void*)&vli[0], (void*)data_, Size*sizeof(BaseInt), cudaMemcpyDeviceToHost ), __LINE__); 					
        }

		operator vli_cpu<BaseInt, Size>() const
        {
            vli_cpu<BaseInt, Size> r;
            copy_vli_to_cpu(r);
            return r;
        }
		
		inline const BaseInt* p() const
        {
            // TODO hide the pointer?
		    return data_;
        }

        inline BaseInt* p()
        {
            // TODO hide the pointer?
            return data_;
        }

		/**
		 multiply and addition operators
		 */
		vli_gpu& operator += (vli_gpu const& vli)
        {
            using vli::plus_assign;
            plus_assign(*this,vli);
            return *this;
        }

		vli_gpu& operator *= (vli_gpu const& vli)
        {
            using vli::multiplies_assign;
            multiplies_assign(*this,vli);
            return *this;
        }

        bool operator == (vli_gpu const& vli) const
        {
            // TODO try on gpu, debug my kernel 
            return vli_cpu<BaseInt, Size>(*this) == vli_cpu<BaseInt, Size>(vli);
        }
        
        bool operator == (vli_cpu<value_type, Size> const& vli) const
        {
             return vli_cpu<value_type, Size> (*this) == vli;
        }
                
        proxy operator[](size_type i) 
        {
            return proxy(data_, i);
        }

        void print(std::ostream& os) const
        {
            os<<vli_cpu<BaseInt, Size>(*this);
        }
		
        void realloc(size_type size) // size = the number of entry of VLI, polynomial or vector
        {
            BaseInt* temp;
            gpu::cu_check_error(cudaMalloc((void**)(temp),size*sizeof(BaseInt)), __LINE__);
            gpu::cu_check_error(cudaMemset((void*)(temp),0, size*sizeof(BaseInt)), __LINE__);	
            gpu::cu_check_error(cudaMemcpy((void*)temp,data_, Size*sizeof(BaseInt),cudaMemcpyDeviceToDevice), __LINE__); 
            std::swap(temp,data_);
            gpu::cu_check_error(cudaFree(data_), __LINE__);
        }
        
        void realloc_lost(size_type size) // datas lost
        {
            BaseInt* temp;
            gpu::cu_check_error(cudaMalloc((void**)(temp),size*sizeof(BaseInt)), __LINE__);
            gpu::cu_check_error(cudaMemset((void*)(temp),0, size*sizeof(BaseInt)), __LINE__);	
            std::swap(temp,data_);
            gpu::cu_check_error(cudaFree(data_), __LINE__);
        }
    
	private:
		BaseInt* data_;
	};
	
	/**
	 multiply and addition operators, suite ...
	 */
	template <class BaseInt, int Size>
	const vli_gpu<BaseInt, Size> operator + (vli_gpu<BaseInt, Size> vli_a, vli_gpu<BaseInt, Size> const& vli_b)
    {
        // TODO check whether direct implementation (without += ) is faster
        vli_a += vli_b;
        return vli_a;
    }
	
	template <class BaseInt, int Size>
	const vli_gpu<BaseInt, Size> operator * (vli_gpu<BaseInt, Size> vli_a, vli_gpu<BaseInt, Size> const& vli_b)
    {
        // TODO check whether direct implementation (without *= ) is faster
        vli_a *= vli_b;
        return vli_a;
    }
	
	template <class BaseInt, int Size>
	std::ostream& operator<<(std::ostream& os, vli_gpu<BaseInt, Size> const& vli)
    {
        vli.print(os);
        return os;
    }
	
}
#endif
