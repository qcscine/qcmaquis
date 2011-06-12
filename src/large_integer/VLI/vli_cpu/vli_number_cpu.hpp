/*
 *  vli.h
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef VLI_NUMBER_CPU_HPP
#define VLI_NUMBER_CPU_HPP
#include <ostream>
#include <vector>
#include "detail/vli_number_cpu_function_hooks.hpp"
#include <cmath>


namespace vli
{
    /**
    template forward declaration 
    */
    template<class BaseInt>
    class vli_gpu;
	
	template<class BaseInt>
    class vli_cpu 
    {
    public:
		typedef BaseInt                       value_type;       // The type T of the elements of the matrix
        typedef BaseInt&                      reference;        // Reference to value_type
        typedef BaseInt const&                const_reference;  // Const reference to value_type
        typedef std::size_t                   size_type;        // Unsigned integer type that represents the dimensions of the matrix
        typedef std::ptrdiff_t                difference_type;  // Signed integer type to represent the distance of two elements in the memory
		
        typedef BaseInt base_int;
        enum { size = SIZE_BITS/(8*sizeof(BaseInt)) };
        
        explicit vli_cpu()
        {
			data_ = (BaseInt*)malloc(size*sizeof(BaseInt));
			memset((void*)data_,0,size*sizeof(BaseInt));
        }
        
		vli_cpu(BaseInt num)
        {
			data_  = (BaseInt*)malloc(size*sizeof(BaseInt));
			memset((void*)data_,0,size*sizeof(BaseInt));			
            *data_ = num;
        }
		
        vli_cpu(vli_cpu const& r)
        {
			data_  = (BaseInt*)malloc(size*sizeof(BaseInt));
            memcpy((void*)data_,(void*)r.data_,size*sizeof(BaseInt));
        }
		
        ~vli_cpu()
        {
			free(data_);
        }
        
        friend void swap(vli_cpu& a, vli_cpu& b)
        {
			std::swap(a.data_,b.data_);
			/*
			*a.data_ ^= *b.data_;
			*b.data_ ^= *a.data_;
			*a.data_ ^= *b.data_;*/			
        }
		
        vli_cpu& operator= (vli_cpu  r)
        {
            swap(*this,r);
            return *this;
        }
        
        BaseInt& operator[](size_type i)
        {
            return *(data_+i);
        }
		
        BaseInt const& operator[](size_type i) const
        {
            return *(data_+i);
        }
		
        /**
		 multiply and addition operators
		 */
        vli_cpu& operator += (vli_cpu const& vli)
        {
            using vli::detail::plus_assign;
            detail::plus_assign(*this,vli);
            return *this;
        }
		
        vli_cpu& operator *= (vli_cpu const& vli)
        {
            using vli::detail::multiplies_assign;
            multiplies_assign(*this,vli);
            return *this;
        }
		
        bool operator == (vli_cpu const& vli) const
        {
			int n = memcmp((void*)data_,(void*)vli.data_,size*sizeof(BaseInt));
			return (0 == n);
        }
        
        void print(std::ostream& os) const
        {
            int i = size - 1 ;
			while( i != 0)
			{
			   	i--;
				os << *(data_+i);
			}
		}
		
		std::size_t BaseTen()
		{
			std::size_t Res = 0;
			for(int i=0;i < size;i++)
				Res+=*(data_+i)*(pow (BASE,i));
			
			return Res;
		}
		
    private:
		BaseInt* data_;
    };
	
    /**
     multiply and addition operators, suite ...
     */
    template <class BaseInt>
    const vli_cpu<BaseInt> operator + (vli_cpu<BaseInt> vli_a, vli_cpu<BaseInt> const& vli_b)
    {
        vli_a += vli_b;
        return vli_a;
    }
    
    template <class BaseInt>
    const vli_cpu<BaseInt> operator * (vli_cpu<BaseInt> vli_a, vli_cpu<BaseInt> const& vli_b)
    {
        vli_a *= vli_b;
        return vli_a;
    }
    
    
    /**
    stream 
    */
    template<typename BaseInt>
    std::ostream& operator<< (std::ostream& os,  vli_cpu<BaseInt> const& vli)
    {
        vli.print(os);
        return os;
    }
     
}

#endif //VLI_NUMBER_CPU_HPP
