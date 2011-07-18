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
#include <ostream>
#include "boost/swap.hpp"
#include "function_hooks/vli_number_cpu_function_hooks.hpp"


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
                data_[i] = r.data_[i];
        }
        
        friend void swap(vli_cpu& a, vli_cpu& b)
        {
			boost::swap(a.data_,b.data_); // gcc 4.2, std::swap does not work with static array
        }
		
        vli_cpu& operator= (vli_cpu r)
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

        void print_raw(std::ostream& os) const
        {
            int i = Size;
            os << "(" ;
			while( i != 0){
			   	i--;
				os << *(data_+i);
                (i == 0) ? (os << ")"):(os << " ");
			}
		}

        void print(std::ostream& os) const
        {
            os<<get_str();
        }

        /**
          * Returns a string with a base10 represenation of the VLI
          */
          
        std::string get_str() const
        {
            vli_cpu tmp(*this);
            if(tmp.is_negative())
            {
                tmp.negate();
                size_type ten_exp = order_of_magnitude_base10(tmp);
                return std::string("-")+get_str_helper_inplace(tmp,ten_exp);
            }
            else
            {
                size_type ten_exp = order_of_magnitude_base10(tmp);
                return get_str_helper_inplace(tmp,ten_exp);
            }
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

        /**
          * Returns the order of magnitude of 'value' in base10
          */
        size_type order_of_magnitude_base10(vli_cpu const& value) const
        {
            assert(!value.is_negative());
            
            vli_cpu value_cpy(value);
            vli_cpu decimal(10);
            size_type exp = 0;

            // Find correct order (10^exp) 
            while(!value_cpy.is_negative())
            {
                value_cpy=value; // reset
                value_cpy-=decimal;
                decimal *= vli_cpu(10);
                ++exp;
            }
            --exp;
            return exp;
        }

        /**
          * A helper function to generate the base10 representation for get_str().
          */
        std::string get_str_helper_inplace(vli_cpu& value, size_type ten_exp) const
        {
            assert(!value.is_negative());

            // Create a number 10^(exponent-1) sin
            vli_cpu dec(1);
            for(size_type e=0; e < ten_exp; ++e)
                dec *= vli_cpu(10);

            // Find the right digit for 10^ten_exp
            vli_cpu value_cpy(value);
            int digit=1;
            while((!value_cpy.is_negative()) && digit<=10)
            {
                value_cpy = value; // reset
                value_cpy-= vli_cpu(digit)*dec;
                ++digit;
            }
            digit-=2; // we got two to far

            assert(digit >=0);
            assert(digit < 10);

            value-= vli_cpu(digit)*dec;

            if(ten_exp <= 0)
                return boost::lexical_cast<std::string>(digit);
            else
                return boost::lexical_cast<std::string>(digit)+get_str_helper_inplace(value,ten_exp-1);
        }


        /**
          * Data members
          */

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
