/*
 *  converter.h
 *
 *  Created by Tim Ewart (timothee.ewart@unige.ch) and  Andreas Hehn (hehn@phys.ethz.ch) on 17.07.11.
 *  Copyright 2011 University of Geneva and Eidgenössische Technische Hochschule Züric. All rights reserved.
 *
 */

#ifndef VLI_CONVERTER_H
#define VLI_CONVERTER_H

#include "boost/lexical_cast.hpp"
namespace vli
{ 

    template<class BaseInt, int Size>
    class converter
    {
    public:
    converter(vli_cpu<BaseInt,Size> const & vli)
    :exponent_old_(-1)
    {
        data_ = new vli_cpu<BaseInt,Size>(vli);
    }
  
    ~converter(){
        delete data_; // think shared pointer for automatic desalocate
    }
  
    const std::string get_string() 
    {
        vli_cpu<BaseInt, Size> tmp(*data_);
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
    
    const std::string  get_string_helper(vli_cpu<BaseInt, Size>& value)
    {
          std::string result;
            // Find correct order (=exponent in (10^exponent) )
            vli_cpu<BaseInt,Size> decimal(1);
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
                decimal *= vli_cpu<BaseInt,Size>(10);
                ++exponent;
            }

            // TODO check
            vli_cpu<BaseInt,Size> dec(1);
            int i;
            
            if(exponent > 0){
                for(size_type e=0; e < exponent-1; ++e)
                    dec *= vli_cpu<BaseInt,Size>(10);
                i = 1;
            }else{
                i = 0;
            }
            // Find digit for 10^exponent
            while(1)
            {
                // TODO int * vli
                value-=vli_cpu<BaseInt,Size>(i)*dec;
      
                if(value.is_negative())
                {
                    // TODO int * vli
                    value+=vli_cpu<BaseInt,Size>(i)*dec;
                    --i;
                    // Subtract the found leading order digit
                    // TODO int * vli
                    value-=vli_cpu<BaseInt,Size>(i)*dec;
                    break;
                }
                value+=vli_cpu<BaseInt,Size>(i)*dec;
                ++i;
            }
     
            if(exponent_old_ != -1 && (exponent_old_ - exponent) > 1) // to avoid to miss a 0 inside the number
            {
                for(int i=0;i< (exponent_old_ - exponent - 1);i++)
                    result+="0";
            }

            result += boost::lexical_cast<std::string>(i);

            if(exponent_old_ == -1 && value.is_null()) // multiple of 10
            {
                for(int e=0 ; e < exponent - e; ++e)
                    result += "0";
   
                return result;         
            }
                                                        
            if(value.is_null()){// finish
                return result;
            }

            exponent_old_= exponent;
            
            result += get_string_helper(value);
            return result;
    }
  
    char const* get_char() 
    {
        return get_string().c_str();
    }
  
   
    void print(std::ostream & os) 
    {
        os << get_string() << std::endl; 
    }
  
    private:
       vli_cpu<BaseInt,Size>* data_;
       std::string name_;
       int exponent_old_;
};

    template<class BaseInt, int Size>
    std::ostream& operator <<( std::ostream & os, converter<BaseInt, Size>  & convert) 
    {
        convert.print(os);
        return os;
    }



}

#endif //VLI_CONVERTER_H