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


#include <ostream>
#include <boost/static_assert.hpp>
#include "vli/detail/bit_masks.hpp"
#include "vli/function_hooks/vli_number_gpu_function_hooks.hpp"
#include "vli/utils/gpu_manager.h"
#include "vli/utils/gpu_base.h"

namespace vli
{
    template<class BaseInt, std::size_t Size>
    class vli_cpu;
    
	template<class BaseInt, std::size_t Size>
	class vli_gpu : public  gpu_array<BaseInt, Size>
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
            
            operator BaseInt() const
            {
                BaseInt base;
                gpu::cu_check_error(cudaMemcpy((void*)&base, (void*)(data_+pos), sizeof(BaseInt),cudaMemcpyDeviceToHost),__LINE__);
                return base;
            }
            
            proxy& operator= (BaseInt i)
            {
                gpu::cu_check_error(cudaMemcpy( (void*)(data_+pos),(void*)&i, sizeof(BaseInt), cudaMemcpyHostToDevice ), __LINE__); 					
                return *this;
            }
            
        private:
            BaseInt* data_;    
            size_type pos;
        };


        
		/**
		constructors 
		*/
        vli_gpu()
        {
        }

		explicit vli_gpu(int num)
        {
            assert( static_cast<BaseInt>((num<0) ? -num : num)  < static_cast<BaseInt>(max_value<BaseInt>::value) );
            
            BaseInt tmp = num & data_mask<BaseInt>::value;
            BaseInt sign = 0x01 & (num>>(sizeof(int)*8-1));
            BaseInt bulk[Size-1];
            for(size_type i = 0; i < Size-2; ++i)
                bulk[i] = sign * data_mask<BaseInt>::value;
            bulk[Size-2] = sign * (data_mask<BaseInt>::value | base<BaseInt>::value);
			
            gpu::cu_check_error(cudaMemcpy((void*)this->data_, (void*)&tmp, sizeof(BaseInt), cudaMemcpyHostToDevice),  __LINE__); //The first number is set to the num value
            gpu::cu_check_error(cudaMemcpy((void*)(this->data_+1), (void*)&bulk[0], (Size-1)*sizeof(BaseInt), cudaMemcpyHostToDevice),  __LINE__);
        }

        explicit vli_gpu(BaseInt* p){
    			gpu::cu_check_error(cudaMemcpy((void*)(this->data_), (void*)p, Size*sizeof(BaseInt), cudaMemcpyDeviceToDevice ), __LINE__);
        }

		vli_gpu(vli_gpu const& vli)
        {
            gpu::cu_check_error(cudaMemcpy((void*)(this->data_), vli.data_, Size*sizeof(BaseInt) , cudaMemcpyDeviceToDevice), __LINE__);
        }

		explicit vli_gpu(vli_cpu<BaseInt, Size> const& vli)
        {
			gpu::cu_check_error(cudaMemcpy((void*)(this->data_), (void*)&vli[0], Size*sizeof(BaseInt), cudaMemcpyHostToDevice ), __LINE__); 		
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
 			gpu::cu_check_error(cudaMemcpy( (void*)&vli[0], (void*)(this->data_), Size*sizeof(BaseInt), cudaMemcpyDeviceToHost ), __LINE__); 					
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
        
        vli_gpu& operator += (BaseInt a)
        {
            using vli::plus_assign;
//          plus_assign(*this,a); doest not work if a is negative to do !
            plus_assign(*this,vli_gpu(a));
            return *this;
        }
        
        vli_gpu& operator -= (vli_gpu const& vli)
        {
            using vli::minus_assign;
            minus_assign(*this,vli);
            return *this;
        }
        
        vli_gpu& operator -= (int a)
        {
            using vli::minus_assign;
            minus_assign(*this,a);
            return *this;
        }

		vli_gpu& operator *= (vli_gpu const& vli)
        {
            // TO DO find a better soluce not fin untill now !
            using vli::multiplies_assign;
            vli_gpu tmp;            
            multiplies_assign(*this,vli,tmp);            
            *this = tmp;
            return *this;
        }
        
        vli_gpu& operator *= (int a)
        {
            using vli::multiplies_assign;
            multiplies_assign(*this, static_cast<BaseInt>(a));
            return *this;
        }
        
        vli_gpu operator - () const
        {
            vli_gpu tmp(*this);
            tmp.negate();
            return tmp;
        }
        
        bool operator == (vli_gpu const& vli) const
        {
            // TODO try on gpu, debug my kernel 
            return (vli_cpu_from(*this) == vli_cpu_from(vli));
        }
        
        bool operator == (vli_cpu<value_type, Size> const& vli) const
        {
             return vli_cpu_from(*this) == vli;
        }

        bool operator < (vli_gpu const& vli) const
        {
            vli_gpu tmp(*this);
            return ( (tmp-=vli).is_negative() );
        }
        
        bool operator < (int i) const
        {
            //TODO improve
            vli_gpu tmp(*this);
            return ( (tmp-=i).is_negative() );
        }
      
        bool operator > (vli_gpu vli) const
        {
            //TODO improve
            return ( (vli-=*this).is_negative() );
        }

        bool operator > (int i) const
        {
            //TODO improve
            vli_gpu tmp(i);
            return ( (tmp-=*this).is_negative() );
        }

        void negate()
        {
            //TODO write test
            using detail::negate;
            using detail::vli_size_tag;
            negate(vli_size_tag<Size>(), this->data_);
        }

        bool is_negative() const
        {
            //TODO write test
            return static_cast<bool>(this->operator[](Size-1)>>data_bits<BaseInt>::value);
        }

        proxy operator[](size_type i) 
        {
            return proxy(this->data_, i);
        }

        const BaseInt operator[] (size_type i) const
        {
            BaseInt base;
            gpu::cu_check_error(cudaMemcpy((void*)&base, (void*)(this->data_+i), sizeof(BaseInt),cudaMemcpyDeviceToHost),__LINE__);
            return base;
        }

        void print(std::ostream& os) const
        {
            os<<vli_cpu_from(*this);
        }

        std::string get_str() const
        {
            return vli_cpu_from(*this).get_str();        
        }
        
    };
	
	/**
	 multiply and addition operators, suite ...
	 */
	template <class BaseInt, std::size_t Size>
	const vli_gpu<BaseInt, Size> operator + (vli_gpu<BaseInt, Size> vli_a, vli_gpu<BaseInt, Size> const& vli_b)
    {
        // TODO check whether direct implementation (without += ) is faster
        vli_a += vli_b;
        return vli_a;
    }
	
    template <class BaseInt, std::size_t Size>
	const vli_gpu<BaseInt, Size> operator + (vli_gpu<BaseInt, Size> vli_a, int b)
    {
        // TODO check whether direct implementation (without += ) is faster
        vli_a += b;
        return vli_a;
    }
    
    template <class BaseInt, std::size_t Size>
	const vli_gpu<BaseInt, Size> operator + (int b, vli_gpu<BaseInt, Size> const& vli_a)
    {
        return vli_a+b;
    }

    template <class BaseInt, std::size_t Size>
	const vli_gpu<BaseInt, Size> operator - (vli_gpu<BaseInt, Size> vli_a, vli_gpu<BaseInt, Size> const& vli_b)
    {
        // TODO check whether direct implementation (without += ) is faster
        vli_a -= vli_b;
        return vli_a;
    }
    
    template <class BaseInt, std::size_t Size>
	const vli_gpu<BaseInt, Size> operator - (vli_gpu<BaseInt, Size> vli_a, int b)
    {
        // TODO check whether direct implementation (without += ) is faster
        vli_a -= b;
        return vli_a;
    }

	template <class BaseInt, std::size_t Size>
	const vli_gpu<BaseInt, Size> operator * (vli_gpu<BaseInt, Size> vli_a, vli_gpu<BaseInt, Size> const& vli_b)
    {
        using vli::multiplies_assign;
        vli_gpu<BaseInt, Size> vli_c;
        multiplies_assign(vli_a, vli_b , vli_c);
        return vli_c;
    }
	
    template <class BaseInt, std::size_t Size>
	const vli_gpu<BaseInt, Size> operator * (vli_gpu<BaseInt, Size> vli_a, int b)
    {        
        vli_a *= b;
        return vli_a;
    }
    
    template <class BaseInt, std::size_t Size>
	const vli_gpu<BaseInt, Size> operator * (int b, vli_gpu<BaseInt, Size> const& vli_a)
    {
        return vli_a*b;
    }
	
	template <class BaseInt, std::size_t Size>
	std::ostream& operator<<(std::ostream& os, vli_gpu<BaseInt, Size> const& vli)
    {
        vli.print(os);
        return os;
    }
	
    /**
      * Ugly conversion functions
      * (due to lack of "explicit" conversion operators in the current standard (will come in C++11))
      */
    template<class BaseInt, std::size_t Size>
    vli_cpu<BaseInt,Size> vli_cpu_from(vli_gpu<BaseInt,Size> const& vli_g)
    {
        vli_cpu<BaseInt, Size> r;
        vli_g.copy_vli_to_cpu(r);
        return r;
    }
    
    template<class BaseInt, std::size_t Size>
    vli_gpu<BaseInt,Size> vli_gpu_from(vli_cpu<BaseInt,Size> const& vli_c)
    {
        return vli_gpu<BaseInt,Size>(vli_c);
    }

}
#endif
