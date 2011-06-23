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

template <class VliType>
class vli_vector_cpu : public std::vector<VliType>
{
public:

	explicit vli_vector_cpu(std::size_t size = 0): std::vector<VliType>(size) // well, std::vector or contiguous mem like GPU ??????
    {
    }
    
	vli_vector_cpu()
    {
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
	vli_vector_cpu& operator *= (VliType const& v)
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
    void print(std::ostream& o) const
    {
        o<<"(";
        typename std::vector<VliType>::const_iterator it = std::vector<VliType>::begin();
        for(; it != std::vector<VliType>::end()-1; ++it)
            o<<*it<<", ";
        o<<*it;
        o<<")";
    }
};
	
/**
 vector addition
 */
template <class VliType>
const vli_vector_cpu<VliType> operator+(vli_vector_cpu<VliType> v_a, vli_vector_cpu<VliType> const& v_b)
{
    v_a += v_b;
    return v_a;
}

/**
  scaling of the vector by an very long integer
  */
template <class VliType>
const vli_vector_cpu<VliType> operator * (vli_vector_cpu<VliType> v_a, VliType const& vli)
{
    v_a *= vli;
    return v_a;
}

template <class VliType>
const vli_vector_cpu<VliType> operator * (VliType const& vli, vli_vector_cpu<VliType> const& v_a)
{
    return v_a * vli;
}

/**
 inner product
 */
template <class VliType>
const VliType inner_product(vli_vector_cpu<VliType> const& v_a, vli_vector_cpu<VliType> const& v_b)
{
    // TODO implement
    assert(v_a.size() == v_b.size());
    return std::inner_product(v_a.begin(), v_a.end(), v_b.begin(), VliType(0));
}


/**
  entrywise multiplication
  */
template <class VliType>
const vli_vector_cpu<VliType> entrywise_product(vli_vector_cpu<VliType> v_a, vli_vector_cpu<VliType> const& v_b)
{
    // TODO implement
    assert(v_a.size() == v_b.size());
    v_a *= entrywise(v_b);
    return v_a;
}

/**
 stream 
 */
template<typename VliType>
std::ostream& operator<< (std::ostream& os,  vli_vector_cpu<VliType> const& v)
{
    v.print(os);
    return os;
}

} //namespace vli

#endif //VLI_VECTOR_CPU_HPP
