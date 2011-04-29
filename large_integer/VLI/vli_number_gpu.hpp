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
#include <stdexcept>
//#include <assert.h>

#include <boost/static_assert.hpp>

#include "boost/lexical_cast.hpp"


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
            //TODO isn't there a better way to initalize the memory on the GPU?
            BaseInt dummy[size];
            for(std::size_t i = 0; i < size; ++i)
                dummy[i] = 0;
		    gpu::check_error( cublasAlloc( size, sizeof(BaseInt), (void**)&data_ ), __LINE__);
            gpu::check_error( cublasSetVector( size, sizeof(BaseInt),&dummy, 1, p(), 1), __LINE__);
        }

		explicit vli_gpu(BaseInt num)
        {
            //TODO isn't there a better way to initalize the memory on the GPU?
            BaseInt dummy[size];
            for(std::size_t i = 0; i < size; ++i)
                dummy[i] = 0;
            gpu::check_error( cublasAlloc( size, sizeof(BaseInt), (void**)&data_ ), __LINE__);
            gpu::check_error( cublasSetVector( size, sizeof(BaseInt), &dummy, 1, p(), 1), __LINE__);
            gpu::check_error( cublasSetVector(1, sizeof(BaseInt), &num, 1, p(), 1), __LINE__);
        }

		vli_gpu(vli_gpu const& vli)
        {
            gpu::check_error( cublasAlloc(size, sizeof(BaseInt), (void**) &data_ ), __LINE__);
            gpu::check_error( cudaMemcpy( data_, vli.data_, size*sizeof(BaseInt) , cudaMemcpyDeviceToDevice), __LINE__);
            std::cout<<"TEST: "<<vli;
            std::cout<<*this;
        }

		vli_gpu(vli_cpu<BaseInt> const& vli)
        {
            // TODO check this!
            BOOST_STATIC_ASSERT( vli_cpu<BaseInt>::size == static_cast<std::size_t>(size) );
		    gpu::check_error(cublasAlloc(size, sizeof(BaseInt), (void**) &data_ ), __LINE__);
		    gpu::check_error(cublasSetVector(vli_cpu<BaseInt>::size, sizeof(BaseInt), &vli[0], 1, p(), 1), __LINE__);
        }

		
		/**
		destructors 
		*/
		~vli_gpu()
        {
            cublasFree(data_);
        }

		//think on the copy swap, I think it is totaly useless for gpu
	//	vli_gpu<T> & operator = (vli_gpu<T> const &  vli);
		/**
		due to gpu architecture the overload is weird
		*/
		void operator = (vli_gpu vli)
        {
            // TODO perhaps the copy-swap implementation is not optimal on the GPU
            swap(vli);
            return *this;
        }

		void swap(vli_gpu& vli)
        {
            using std::swap;
            swap(data_, vli.data_);
        }

		void copy_vli_to_cpu(vli::vli_cpu<BaseInt>& vli) const
        {
            BOOST_STATIC_ASSERT( vli_cpu<BaseInt>::size == static_cast<std::size_t>(size) );
            gpu::check_error(cublasGetVector(size, sizeof(BaseInt), data_  ,1, &vli[0], 1), __LINE__);			
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

        void print(std::ostream& os) const
        {
            os<<vli_cpu<BaseInt>(*this);
        }
		
	private:
		base_int_type* data_;
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
        return vli_b;
    }
	
	template <class BaseInt>
	std::ostream& operator<<(std::ostream& os, vli_gpu<BaseInt> const& vli)
    {
        vli.print(os);
        return os;
    }
	
}
#endif
