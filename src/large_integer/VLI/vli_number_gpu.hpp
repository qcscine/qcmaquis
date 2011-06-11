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
#include "detail/vli_number_gpu_function_hooks.hpp"
#include "GpuManager.h"

namespace vli
{

		
	template<class BaseInt>
	class vli_gpu 
	{
	public:
        typedef BaseInt base_int_type;
        typedef std::size_t size_type;
        enum { size = SIZE_BITS/(8*sizeof(BaseInt)) };

		/**
		constructors 
		*/
		vli_gpu()
        {
		    gpu::cu_check_error(cudaMalloc((void**)&data_, size*sizeof(BaseInt)), __LINE__);			
	    	gpu::cu_check_error(cudaMemset((void*)data_, 0, size*sizeof(BaseInt)), __LINE__);			
        }

		explicit vli_gpu(BaseInt num)
        {
            gpu::cu_check_error(cudaMalloc((void**)&data_, size*sizeof(BaseInt)), __LINE__);			
			gpu::cu_check_error(cudaMemset((void*)data_, 0, size*sizeof(BaseInt)), __LINE__); //Everything is set to 0	
			gpu::cu_check_error(cudaMemset((void*)data_, num, sizeof(BaseInt)), __LINE__); //The first number is set to the num value			
        }

		vli_gpu(vli_gpu const& vli)
        {
			gpu::cu_check_error(cudaMalloc((void**)&data_, size*sizeof(BaseInt)), __LINE__);
            gpu::cu_check_error(cudaMemcpy((void*)data_, vli.data_, size*sizeof(BaseInt) , cudaMemcpyDeviceToDevice), __LINE__);
        }

		vli_gpu(vli_cpu<BaseInt> const& vli)
        {
            BOOST_STATIC_ASSERT(vli_cpu<BaseInt>::size == static_cast<std::size_t>(size) );
			gpu::cu_check_error(cudaMalloc((void**)&data_, size*sizeof(BaseInt)  ), __LINE__);
			gpu::cu_check_error(cudaMemcpy((void*)data_, (void*)&vli[0], vli_cpu<BaseInt>::size*sizeof(BaseInt), cudaMemcpyHostToDevice ), __LINE__); 		
        }
		
		/**
		destructors 
		*/
		~vli_gpu()
        {
            cublasFree(data_);
        }

		friend void swap(vli_gpu& a, vli_gpu& b)
        {
            using std::swap;
            swap(a.data_, b.data_);
        }

		vli_gpu& operator = (vli_gpu vli)
        {
            // TODO perhaps the copy-swap implementation is not optimal on the GPU
            swap(*this,vli);
            return *this;
        }

		void copy_vli_to_cpu(vli::vli_cpu<BaseInt>& vli) const
        {
            BOOST_STATIC_ASSERT( vli_cpu<BaseInt>::size == static_cast<std::size_t>(size) );
 			gpu::cu_check_error(cudaMemcpy( (void*)&vli[0], (void*)data_, vli_cpu<BaseInt>::size*sizeof(BaseInt), cudaMemcpyDeviceToHost ), __LINE__); 					
        }

		operator vli_cpu<BaseInt>() const
        {
            vli_cpu<BaseInt> r;
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
            using vli::detail::plus_assign;
            plus_assign(*this,vli);
            return *this;
        }

		vli_gpu& operator *= (vli_gpu const& vli)
        {
            using vli::detail::multiplies_assign;
            multiplies_assign(*this,vli);
            return *this;
        }

        bool operator == (vli_gpu const& vli) const
        {
            // TODO this should also work directly on the gpu
            return vli_cpu<BaseInt>(*this) == vli_cpu<BaseInt>(vli);
        }

        void print(std::ostream& os) const
        {
            os<<vli_cpu<BaseInt>(*this);
        }
		
	private:
		BaseInt* data_;
	};
	
	/**
	 multiply and addition operators, suite ...
	 */
	template <class BaseInt>
	const vli_gpu<BaseInt> operator + (vli_gpu<BaseInt> vli_a, vli_gpu<BaseInt> const& vli_b)
    {
        // TODO check whether direct implementation (without += ) is faster
        vli_a += vli_b;
        return vli_a;
    }
	
	template <class BaseInt>
	const vli_gpu<BaseInt> operator * (vli_gpu<BaseInt> vli_a, vli_gpu<BaseInt> const& vli_b)
    {
        // TODO check whether direct implementation (without *= ) is faster
        vli_a *= vli_b;
        return vli_a;
    }
	
	template <class BaseInt>
	std::ostream& operator<<(std::ostream& os, vli_gpu<BaseInt> const& vli)
    {
        vli.print(os);
        return os;
    }
	
}
#endif
