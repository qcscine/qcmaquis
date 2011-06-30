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
#include "boost/swap.hpp"
//#include "engine/engine.h"

namespace vli
{
    /**
    template forward declaration 
    */    
    template<class BaseInt, int Size>
    class vli_gpu;
	
	template<class BaseInt, int Size>
    class vli_cpu 
    {
    public:
	    typedef BaseInt                       value_type;       // The type T of the elements of the matrix
        typedef BaseInt&                      reference;        // Reference to value_type
        typedef BaseInt const&                const_reference;  // Const reference to value_type
        typedef std::size_t                   size_type;        // Unsigned integer type that represents the dimensions of the matrix
        typedef std::ptrdiff_t                difference_type;  // Signed integer type to represent the distance of two elements in the memory
        typedef BaseInt base_int;
        
        vli_cpu()
        {
            for(int i=0; i<Size; ++i)
                data_[i] = 0;
        }
        
		explicit vli_cpu(BaseInt num)
        {
            data_[0] = num;
            for(int i=1; i<Size; ++i)
                data_[i]=0;
        }
		
        vli_cpu(vli_cpu const& r)
        {
            for(int i=0; i<Size; ++i)
                data_[i] = r.data_[i]; // because copy constructor
        }
        
        vli_cpu(BaseInt* p)
        {
            for(int i=0; i<Size; ++i)
                data_[i] = p->data_[i]; // because copy constructor
        }
        
        friend void swap(vli_cpu& a, vli_cpu& b)
        {
			boost::swap(a.data_,b.data_); // gcc 4.2, std::swap does not work with static array
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
		
        const BaseInt& operator[](size_type i) const
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
            detail:multiplies_assign(*this,vli);
            return *this;
        }
		
        bool operator == (vli_cpu const& vli) const
        {
			int n = memcmp((void*)data_,(void*)vli.data_,Size*sizeof(BaseInt));
			return (0 == n);
        }
        
        void print(std::ostream& os) const
        {
            int i = Size;
            os << "(" ;
			while( i != 0){
			   	i--;
				os << *(data_+i);
                (i == 0) ? (os << ")"):(os << " ");
			}
		}
		
		std::size_t BaseTen()
		{
			std::size_t Res = 0;
			for(int i=0;i < Size;i++)
				Res+=data_[i]*(pow (BASE,i));
			
			return Res;
		}
		
    private:
		BaseInt data_[Size];
    };
	
    /**
     multiply and addition operators, suite ...
     */
    template <class BaseInt, int Size>
    const vli_cpu<BaseInt, Size> operator + (vli_cpu<BaseInt, Size> vli_a, vli_cpu<BaseInt, Size> const& vli_b)
    {
        vli_a += vli_b;
        return vli_a;
    }
    
    template <class BaseInt, int Size>
    const vli_cpu<BaseInt, Size> operator * (vli_cpu<BaseInt, Size> vli_a, vli_cpu<BaseInt, Size> const& vli_b)
    {
        vli_a *= vli_b;
        return vli_a;
    }
    
    /**
    stream 
    */
    template<typename BaseInt, int Size>
    std::ostream& operator<< (std::ostream& os,  vli_cpu<BaseInt, Size> const& vli)
    {
        vli.print(os);
        return os;
    }
     
}

#endif //VLI_NUMBER_CPU_HPP
