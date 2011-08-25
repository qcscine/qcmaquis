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
#include <detail/bit_masks.hpp>
#include <boost/lexical_cast.hpp>
#include <ostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstring>
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
        
        explicit vli_cpu(int num)
        {            
            assert( static_cast<BaseInt>((num<0) ? -num : num)  < static_cast<BaseInt>(max_value<BaseInt>::value) );
            data_[0] = num & data_mask<BaseInt>::value;
            BaseInt sign = 0x01 & (num>>(sizeof(int)*8-1));
            for(int i=1; i<Size-1; ++i)
                data_[i] = sign * data_mask<BaseInt>::value;
            data_[Size-1] = sign * (base<BaseInt>::value+data_mask<BaseInt>::value);
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
        
        
        vli_cpu& operator += (vli_cpu const& vli)
        {
            using vli::plus_assign;
            plus_assign(*this,vli);
            return *this;
        }
        
        vli_cpu& operator += (int a)
        {
            using vli::plus_assign;
            //  plus_assign(*this,a); doest not work if a is negative to do ! 
            plus_assign(*this,vli_cpu(a));
            return *this;
        }
        
        vli_cpu& operator -= (vli_cpu const& vli)
        {
            using vli::plus_assign;
            vli_cpu tmp(vli);
            tmp.negate();
            plus_assign(*this,tmp);
            return *this;
        }

        vli_cpu& operator -= (int a)
        {
            using vli::plus_assign;
            plus_assign(*this,vli_cpu(-a));
            return *this;
        }

        vli_cpu& operator *= (int a)
        {
            assert(true);// add a test if a is very large
            using vli::multiplies_assign;
    
            int na(1),nb(1);        
                  
            if((*this).is_negative()){
                (*this).negate();
                na = -1;
            }
        
            if(a<0){
                a *=-1;
                nb = -1;
            }
        
            multiplies_assign(*this,a);
        
            if(nb*na == -1)
                (*this).negate();
       
            return *this;
        }

        vli_cpu& operator *= (vli_cpu const& a)
        {
            using vli::multiplies_assign;
            bool result_is_negative = static_cast<bool>((this->data_[Size-1] ^ a[Size-1]) >> data_bits<BaseInt>::value);
            if(result_is_negative)// test if 
            {
                this->negate(); // - to +
                multiplies_assign(*this,a); // +*- or -*+ 
                this->negate(); // + to -
            }
            else
            {
                multiplies_assign(*this,a); // +*+ or -*- the default multiplication works
            }
            return *this;
        }

        vli_cpu operator - () const
        {
            vli_cpu tmp(*this);
            tmp.negate();
            return tmp;
        }

        bool operator == (vli_cpu const& vli) const
        {
			int n = memcmp((void*)data_,(void*)vli.data_,Size*sizeof(BaseInt));
			return (0 == n);
        }
        
        bool operator != (vli_cpu const& vli) const
        {
           for(std::size_t i(0); i < Size-1 ; ++i)
           {
               if((*this)[i] != vli[i])
                   return true;
           }        
           return false;
        }
        

        bool operator < (vli_cpu const& vli) const
        {
            // TODO improve
            vli_cpu tmp(*this);
            return ( (tmp-=vli).is_negative() );
        }

        bool operator < (int i) const
        {
            //TODO improve
            vli_cpu tmp(*this);
            return ( (tmp-=i).is_negative() );
        }
      
        bool operator > (vli_cpu vli) const
        {
            //TODO improve
            return ( (vli-=*this).is_negative() );
        }

        bool operator > (int i) const
        {
            //TODO improve
            vli_cpu tmp(i);
            return ( (tmp-=*this).is_negative() );
        }

        void negate()
        {
            for(size_type i=0; i < Size-1; ++i)
                data_[i] = (~data_[i])&data_mask<BaseInt>::value;
            data_[Size-1] = (~data_[Size-1])&(base<BaseInt>::value+data_mask<BaseInt>::value);
            (*this)+=vli_cpu(1);
        }

        bool is_negative() const
        {
            return static_cast<bool>((data_[Size-1]>>data_bits<BaseInt>::value));
        }
        
        void print_raw(std::ostream& os) const
        {
            os << "(" ;
            for(std::size_t i=Size-1; i > 0; --i)
				os << std::showbase << std::hex << data_[i]<<" ";

		    os << std::showbase << std::hex << data_[0];
            os << ")";
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
            
            vli_cpu<BaseInt, size*2> tmp;
            
            if((*this).is_negative())
            {
                const_cast<vli_cpu<BaseInt, Size> & >(*this).negate();
                
                for(int i=0; i<size; ++i)
                    tmp[i] = (*this)[i];
                 
                tmp.negate();
                const_cast<vli_cpu<BaseInt, Size> & >(*this).negate();
            }else{
                
                for(int i=0; i<size; ++i)
                    tmp[i] = (*this)[i];
                
            }
                    
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

    private:

        /**
          * Returns the order of magnitude of 'value' in base10
          */
        template<int SizeDouble>
        size_type order_of_magnitude_base10(vli_cpu<BaseInt,SizeDouble> const& value) const
        {
            assert(!value.is_negative());
            
            vli_cpu<BaseInt,SizeDouble> value_cpy(value);
            vli_cpu<BaseInt,SizeDouble> decimal(1);
            size_type exp = 0;

            // Find correct order (10^exp) 
            while(!value_cpy.is_negative())
            {
                value_cpy=value; // reset
                vli_cpu<BaseInt,SizeDouble> previous_decimal(decimal);
                decimal *= 10;
                ++exp;
                if(decimal < previous_decimal) // Overflow! (we can handle it.)
                {
                    break;
                }
                value_cpy-=decimal;
            }
            --exp;
            return exp;
        }

        /**
          * A helper function to generate the base10 representation for get_str().
          */
        template<int SizeDouble>
        std::string get_str_helper_inplace(vli_cpu<BaseInt,SizeDouble>& value, size_type ten_exp) const
        {
            assert(!value.is_negative());

            // Create a number 10^(exponent-1) sin
            vli_cpu<BaseInt,SizeDouble> dec(1);
            for(size_type e=0; e < ten_exp; ++e)
                dec *= 10;

            // Find the right digit for 10^ten_exp
            vli_cpu<BaseInt,SizeDouble> value_cpy(value);
            int digit=0;
            while((!value_cpy.is_negative()) && digit<=11)
            {
                value_cpy = value; // reset
                ++digit;
                if(digit*dec < (digit-1)*dec) // Overflow (we can handle it.)
                {
                    break;
                }
                value_cpy-= digit*dec;
            }
            --digit; // we went to far

            assert(digit >=0);
            assert(digit < 10); 

            value-= digit*dec;

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
    const vli_cpu<BaseInt, Size> operator + (vli_cpu<BaseInt, Size> vli_a, int b)
    {
        vli_a += b;
        return vli_a;
    }
    
    template <class BaseInt, int Size>
    const vli_cpu<BaseInt, Size> operator + (int b, vli_cpu<BaseInt, Size> const& vli_a)
    {
        return vli_a+b;
    }
    
    template <class BaseInt, int Size>
    const vli_cpu<BaseInt, Size> operator - (vli_cpu<BaseInt, Size> vli_a, vli_cpu<BaseInt, Size> const& vli_b)
    {
        vli_a -= vli_b;
        return vli_a;
    }
    
    template <class BaseInt, int Size>
    const vli_cpu<BaseInt, Size> operator - (vli_cpu<BaseInt, Size> vli_a, int b)
    {
        vli_a -= b;
        return vli_a;
    }

    template <class BaseInt, int Size>
    const vli_cpu<BaseInt, Size> operator * (vli_cpu<BaseInt, Size>  vli_a, vli_cpu<BaseInt, Size> const& vli_b)
    {
        vli_a *= vli_b;
        return vli_a;
    /*
        vli_cpu<BaseInt, 2*Size> vli_res;
        
        int na(1),nb(1);        
                
        if(vli_a.is_negative()){
            const_cast<vli_cpu<BaseInt, Size> & >(vli_a).negate();
            na = -1;
        }
        
        if(vli_b.is_negative()){
            const_cast<vli_cpu<BaseInt, Size> & >(vli_b).negate();
            nb = -1;
        }

        multiplies<BaseInt, 2*Size, Size>(vli_res, vli_a, vli_b);
        
        if(nb*na == -1)
            vli_res.negate();
       
        if(na == -1){
           const_cast<vli_cpu<BaseInt, Size> & >(vli_a).negate();   
        }

        if(nb == -1){
           const_cast<vli_cpu<BaseInt, Size> & >(vli_b).negate();   
        }

        return vli_res;
    */
    }

    template <class BaseInt, int Size>
    const vli_cpu<BaseInt, Size> operator * (vli_cpu<BaseInt, Size> vli_a, int b)
    {
        vli_a *= b;
        return vli_a;
    }

    template <class BaseInt, int Size>
    const vli_cpu<BaseInt, Size> operator * (int b, vli_cpu<BaseInt, Size> const& a)
    {
        return a*b;
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
