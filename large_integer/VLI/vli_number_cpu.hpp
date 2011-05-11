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
        typedef BaseInt base_int;
        typedef std::size_t size_type;
        enum { size = SIZE_BITS/(8*sizeof(BaseInt)) };
        
        /**
        constructors 
        */
        vli_cpu()
            : data_(size)
        {
        }
        
        explicit vli_cpu(BaseInt num)
            : data_(size)
        {
            data_[0] = num;
        }

        vli_cpu(vli_cpu const& r)
            :data_(r.data_)
        {
        }

        /**
        destructor 
        */
        ~vli_cpu()
        {
        }
        
        /**
          swap
          */
        friend void swap(vli_cpu& a, vli_cpu& b)
        {
            swap(a.data_,b.data_);
        }

        /**
          assignment operator
          */
        vli_cpu& operator= (vli_cpu  r)
        {
            swap(*this,r);
            return *this;
        }
        
        /**
        logistics operators 
        */
        BaseInt& operator[](size_type i)
        {
            return data_[i];
        }

        BaseInt const& operator[](size_type i) const
        {
            return data_[i];
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
            return (data_ == vli.data_);
        }
        
        void print(std::ostream& os) const
        {
			typename std::vector<BaseInt>::const_iterator it = data_.end();	
            /** I reversed the print easier to read **/
			while (it != data_.begin())
			{
				it--;
				os << *it << " " ;
			}
		}
		
		std::size_t BaseTen()
		{
			std::size_t Res = 0;
			for(int i=0;i < size;i++)
               Res+=data_[i]*(pow (BASE,i));
			
			return Res;
		}
						  
    private:
        std::vector<BaseInt> data_;

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
