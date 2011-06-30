/*
 *  vli_vector_cpu.h
 *  Untitled
 *
 *  Created by Tim Ewart on 23/03/11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#ifndef VLI_VECTOR_CPU_HPP
#define VLI_VECTOR_CPU_HPP

#include "detail/vli_vector_cpu_function_hooks.hpp"
#include "vli_utils/utils.hpp"
#include <numeric>

namespace vli
{

template <class VLI>
class vli_vector_cpu 
{
public:
    
    typedef typename VLI::value_type value_type;
    typedef typename VLI::vli_size vli_size;    
    typedef std::size_t size_type;
    
    explicit vli_vector_cpu(size_type size = 1): size_(size){
        //to have a better compatibility with the gpu vector
        data_ = (BaseInt*)malloc(vli_size*size_*sizeof(BaseInt));
        memset((void*)data_,0,vli_size*size_*sizeof(BaseInt));
    }
    
    vli_vector_cpu(vli_vector_cpu<VLI> const & r):size_(r.size()){
        assert( this->size() == r.size() );
        data_ = (BaseInt*)malloc(vli_size*size_*sizeof(BaseInt));
        memcpy((void*)data_, (void*)r[0], vli_size*size_*sizeof(BaseInt));
    }
    
    friend void swap(vli_vector_cpu  &a, vli_vector_cpu  &b){
        using std::swap;
        swap(a.data_,b.data_);
        swap(a.size_,b.size_);
    }
    
    BaseInt* operator[](size_type i){ //return the position of the VLI
        return (data_+i*vli_size);
    }
    
    void init(size_type i, size_type amount, ...){
        assert(amount <= vli_size);
        va_list vl;
        va_start(vl,amount);
        for (int j=0;j<amount;j++){
            BaseInt val=va_arg(vl,size_type);
            *(data_+i*vli_size+j) = val;
        }
        va_end(vl);
    }
    
    const BaseInt* operator[](size_type i) const {
        return (data_+i*vli_size);
    }
        
    /**
    getting the number of element inside the vector
    */
    size_type size() const {
        return size_;
    }
 
    size_type size(){
        return size_;
    }
    /**
     operator =
     */    
    vli_vector_cpu & operator=(vli_vector_cpu v){
        swap(*this, v);
        return *this;
    }
    
	/**
	 plus assign (vector addition)
	 */															
	vli_vector_cpu& operator += (vli_vector_cpu const& v)
    {
        using detail::plus_assign;
        plus_assign(*this, v);
        return *this;
    }

    /**
      multiply by a scalar
      */
	vli_vector_cpu& operator *= (VLI const& v)
    {
        using detail::multiplies_assign;
        multiplies_assign(*this, v);
        return *this;
    }
    
    /**
      entrywise multiplication
      */
    vli_vector_cpu& operator *= (utils::entrywise<vli_vector_cpu> const& v)
    {
        using detail::entrywise_multiplies_assign;
        entrywise_multiplies_assign(*this, v.vector);
        return *this;
    }

    /**
      print to ostream
      */
    void print(std::ostream& os) const
    {
        os<<"(";
        for(size_type i=0; i < size_ ; i++){
            int j = vli_size - 1 ;
			while( j != 0) {
			   	j--;
				os << *((data_+i*vli_size+j)) << " ";
			}
            i==(size_-1) ? (os << ")") : (os << "|");
        }
        os << std::endl;
    }
    
private:
    BaseInt* data_;
    size_type size_;
};
	
/**
 vector addition
 */
template <class VLI>
const vli_vector_cpu<VLI> operator+(vli_vector_cpu<VLI> v_a, vli_vector_cpu<VLI> const& v_b)
{
    v_a += v_b;
    return v_a;
}

/**
  scaling of the vector by an very long integer
  */
template <class VLI>
const vli_vector_cpu<VLI> operator * (vli_vector_cpu<VLI> v_a, VLI const& vli)
{
    v_a *= vli;
    return v_a;
}

template <class VLI>
const vli_vector_cpu<VLI> operator * (VLI const& vli, vli_vector_cpu<VLI> const& v_a)
{
    return v_a * vli;
}

/**
 inner product
 */
template <class VLI>
const VLI inner_product(vli_vector_cpu<VLI> const& v_a, vli_vector_cpu<VLI> const& v_b)
{
    // TODO implement
    assert(v_a.size() == v_b.size());
    return std::inner_product(v_a.begin(), v_a.end(), v_b.begin(), VLI(0));
}


/**
  entrywise multiplication
  */
template <class VLI>
const vli_vector_cpu<VLI> entrywise_product(vli_vector_cpu<VLI> v_a, vli_vector_cpu<VLI> const& v_b)
{
    // TODO implement
    assert(v_a.size() == v_b.size());
    v_a *= entrywise(v_b);
    return v_a;
}

/**
 stream 
 */
template<typename VLI>
std::ostream& operator<< (std::ostream& os,  vli_vector_cpu<VLI> const& v)
{
    v.print(os);
    return os;
}

} //namespace vli

#endif //VLI_VECTOR_CPU_HPP
