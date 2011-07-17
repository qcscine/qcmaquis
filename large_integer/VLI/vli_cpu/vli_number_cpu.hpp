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
#include <boost/lexical_cast.hpp>
#include <ostream>
#include <vector>
#include <string>
#include <cmath>
#include "boost/swap.hpp"
#include "function_hooks/vli_number_cpu_function_hooks.hpp"
#include <iostream>

    


namespace vli
{
    /**
    template forward declaration 
    */    
	
	template<class BaseInt, int Size>
    class vli_cpu 
    {
    public:
	    typedef BaseInt         value_type;     // Data type to store parts of the very long integer (usually int)
        typedef std::size_t     size_type;      // Size type of the very long integers (number of parts)
        enum {size = Size};                     // Number of parts of the very long integer (eg. how many ints)
        
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
        
//        vli_cpu(BaseInt* p)
//        {
//            for(int i=0; i<Size; ++i)
//                data_[i] = p->data_[i]; // because copy constructor
//        }
        
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
		 multiply,substraction, addition operators
		 */
        
        vli_cpu& operator -= (vli_cpu const& vli)
        {
            using vli::plus_assign;
            vli_cpu tmp(vli);
            tmp.negate();
            plus_assign(*this,tmp);
            return *this;
        }
        
        vli_cpu& operator += (vli_cpu const& vli)
        {
            using vli::plus_assign;
            plus_assign(*this,vli);
            return *this;
        }
		
        vli_cpu& operator *= (vli_cpu const& vli)
        {
            using vli::multiplies_assign;
            multiplies_assign(*this,vli);
            return *this;
        }
		
        bool operator == (vli_cpu const& vli) const
        {
			int n = memcmp((void*)data_,(void*)vli.data_,Size*sizeof(BaseInt));
			return (0 == n);
        }
        
        void negate()
        {
            for(size_type i=0; i < Size-1; ++i)
                data_[i] = (~data_[i])&BASE_MINUS;
           data_[Size-1] = (~data_[Size-1])&(BASE+BASE_MINUS);//0x1FF;
           (*this)+=vli_cpu(1);
        }

        bool is_negative() const
        {
            return static_cast<bool>((data_[Size-1]>>LOG_BASE));
        }
        
        bool is_null() const
        {
            for(int i=0; i < Size - 1 ; ++i){
                if(data_[i] != 0)
                    return false;
            }
            return true;
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

		char const* get_char() const
        {
            return get_string().c_str();
        }
        
        
        std::string get_string() const
        {
            vli_cpu tmp(*this);
            if(tmp.is_negative())
            {
                tmp.negate();
                return std::string("-")+get_string_helper(tmp);
            }
            else
            {
                return get_string_helper(tmp);
            }
        }

        // TODO make private
        std::string get_string_helper(vli_cpu& value) const
        {
            std::string result;
            // Find correct order (=exponent in (10^exponent) )
            vli_cpu decimal(1);
            int exponent = 0;
            while(1)
            {
                value-=decimal;
                if(value.is_negative())
                {
                    value+=decimal;
                    break;
                }
                value+=decimal;
                decimal *= vli_cpu(10);
                ++exponent;
            }

            // TODO check
            vli_cpu dec(1);
            int i;
            
            if(exponent > 0){
                for(size_type e=0; e < exponent-1; ++e)
                    dec *= vli_cpu(10);
                i = 1;
            }else{
                i = 0;
            }
            // Find digit for 10^exponent
            while(1)
            {
                // TODO int * vli
                value-=vli_cpu(i)*dec;
      
                if(value.is_negative())
                {
                    // TODO int * vli
                    value+=vli_cpu(i)*dec;
                    --i;
                    // Subtract the found leading order digit
                    // TODO int * vli
                    value-=vli_cpu(i)*dec;
                    break;
                }
                value+=vli_cpu(i)*dec;
                ++i;
            }
            
            if(exponent <= 1)
            {
                 return boost::lexical_cast<std::string>(i);
            }
            
            result += boost::lexical_cast<std::string>(i);

            if(value.is_null()){// if multiple of 0
                for(int e=0 ; e < exponent - e; ++i)
                    result += boost::lexical_cast<std::string>(0);

                return result;
            }

            result += get_string_helper(value);
            return result;
        }
        
		long int BaseTen() // for debuging on small BASE
		{
			long int Res = 0;
            if(this->is_negative())
            {
                this->negate();
                for(size_type i=0;i < Size;i++)
                    Res+=(data_[i]&BASE_MINUS)*(pow(BASE,i));
                Res *= -1;
                this->negate();
            }
            else
            {
                for(size_type i=0;i < Size;i++)
                    Res+=(data_[i]&BASE_MINUS)*(pow(BASE,i));
            }
            return Res;
		}
        
    private:
		BaseInt data_[Size] __attribute__ ((aligned (16)));
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
