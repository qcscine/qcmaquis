/*
 6*  vli.h
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef VLI_NUMBER_CPU_HPP
#define VLI_NUMBER_CPU_HPP
#include "vli/detail/bit_masks.hpp"
#include "vli/function_hooks/vli_number_cpu_function_hooks.hpp"
#include <boost/lexical_cast.hpp>
#include <vector>
#include <string>
#include <cassert>
#include <cstring>
#include <ostream>
#include <sstream>
#include <boost/swap.hpp>


namespace vli{
    
    template<class BaseInt, std::size_t Size>
    class vli_cpu 
    {
    public:
        typedef BaseInt         value_type;     // Data type to store parts of the very long integer (usually int) -
        typedef std::size_t     size_type;      // Size type of the very long integers (number of parts)
        enum {size = Size};                     // Number of parts of the very long integer (eg. how many ints)
        // c - constructors, copy-swap, access   
        vli_cpu();
        explicit vli_cpu(int num);
        vli_cpu(vli_cpu const& r);
        void swap(vli_cpu& a, vli_cpu& b);
        vli_cpu& operator= (vli_cpu r);
        BaseInt& operator[](size_type i);
        const BaseInt& operator[](size_type i) const;
        // c - negative number
        void negate();
        bool is_negative() const;
        // c - basic operator
        vli_cpu& operator += (vli_cpu const& vli);        
        vli_cpu& operator += (BaseInt const a);
        vli_cpu& operator -= (vli_cpu const& vli);
        vli_cpu& operator -= (BaseInt a);
        vli_cpu& operator *= (BaseInt a); // 192 bits -> 192 bits
        vli_cpu& operator *= (vli_cpu const& a); // 192 bits -> 192 bits
        
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
           for(size_type i(0); i < size-1 ; ++i)
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

        void print_raw(std::ostream& os) const;
        void print(std::ostream& os) const;
        
        std::string get_str() const;
        size_type order_of_magnitude_base10(vli_cpu<BaseInt,size> const& value) const;
        std::string get_str_helper_inplace(vli_cpu<BaseInt,size>& value, size_type ten_exp) const;
    private:
        BaseInt data_[Size];
    };
    
    /**
     multiply and addition operators, suite ...
     */
    template <class BaseInt, std::size_t Size>
    const vli_cpu<BaseInt, Size> operator + (vli_cpu<BaseInt, Size> vli_a, vli_cpu<BaseInt, Size> const& vli_b)
    {
        vli_a += vli_b;
        return vli_a;
    }
    
    template <class BaseInt, std::size_t Size>
    const vli_cpu<BaseInt, Size> operator + (vli_cpu<BaseInt, Size> vli_a, int b)
    {
        vli_a += b;
        return vli_a;
    }
    
    template <class BaseInt, std::size_t Size>
    const vli_cpu<BaseInt, Size> operator + (int b, vli_cpu<BaseInt, Size> const& vli_a)
    {
        return vli_a+b;
    }
    
    template <class BaseInt, std::size_t Size>
    const vli_cpu<BaseInt, Size> operator - (vli_cpu<BaseInt, Size> vli_a, vli_cpu<BaseInt, Size> const& vli_b)
    {
        vli_a -= vli_b;
        return vli_a;
    }
    
    template <class BaseInt, std::size_t Size>
    const vli_cpu<BaseInt, Size> operator - (vli_cpu<BaseInt, Size> vli_a, int b)
    {
        vli_a -= b;
        return vli_a;
    }

    template <class BaseInt, std::size_t Size>
    void multi_nt(vli_cpu<BaseInt, 2*Size>& vli_res, vli_cpu<BaseInt, Size> const&  vli_a, vli_cpu<BaseInt, Size> const& vli_b) // C nt = non truncated
    {         
        int na(1),nb(1);        
         
        if(vli_a.is_negative()){
        const_cast<vli_cpu<BaseInt, Size> & >(vli_a).negate();
        na = -1;
        }
         
        if(vli_b.is_negative()){
        const_cast<vli_cpu<BaseInt, Size> & >(vli_b).negate();
        nb = -1;
        }
         
        multiplies<BaseInt, Size>(vli_res, vli_a, vli_b);
         
        if(nb*na == -1)
        vli_res.negate();
         
        if(na == -1){
        const_cast<vli_cpu<BaseInt, Size> & >(vli_a).negate();   
        }
         
        if(nb == -1){
        const_cast<vli_cpu<BaseInt, Size> & >(vli_b).negate();   
        }
    }
    
    template <class BaseInt, std::size_t Size>
    const vli_cpu<BaseInt, Size> operator * (vli_cpu<BaseInt, Size>  vli_a, vli_cpu<BaseInt, Size> const& vli_b)
    {
        vli_a *= vli_b;
        return vli_a;
    }

    template <class BaseInt, std::size_t Size>
    const vli_cpu<BaseInt, Size> operator * (vli_cpu<BaseInt, Size> vli_a, int b)
    {
        vli_a *= b;
        return vli_a;
    }

    template <class BaseInt, std::size_t Size>
    const vli_cpu<BaseInt, Size> operator * (int b, vli_cpu<BaseInt, Size> const& a)
    {
        return a*b;
    }
    
    /**
    stream 
    */
    template<typename BaseInt, std::size_t Size>
    std::ostream& operator<< (std::ostream& os,  vli_cpu<BaseInt, Size> const& vli);
  
#include "vli/vli_cpu.hpp"
    
}

#endif //VLI_NUMBER_CPU_HPP
