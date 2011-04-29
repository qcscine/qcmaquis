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

namespace vli
{

template <class VliType>
class vli_vector : public std::vector<VliType>
{
public:

	vli_vector(std::size_t size)
        : std::vector<VliType>(size)
    {
    }
	~vli_vector()
    {
    }

	/**
	 multiply and addition operators
	 */															
	vli_vector& operator += (vli_vector const& v)
    {
        using detail::plus_assign;
        plus_assign(*this, v);
        return *this;
    }

	vli_vector& operator *= (vli_vector const& v)
    {
        using detail::multiplies_assign;
        multiplies_assign(*this, v);
        return *this;
    }

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
 multiply and addition operators, suite ...
 */
template <class VliType>
const vli_vector<VliType> operator+(vli_vector<VliType> v_a, vli_vector<VliType> const& v_b)
{
    v_a += v_b;
    return v_a;
}

template <class VliType>
const vli_vector<VliType> operator*(vli_vector<VliType> v_a, vli_vector<VliType> const& v_b)
{
    v_a *= v_b;
    return v_a;
}


/**
 stream 
 */
template<typename VliType>
std::ostream& operator<< (std::ostream& os,  vli_vector<VliType> const& v)
{
    v.print(os);
    return os;
}

} //namespace vli

#endif //VLI_VECTOR_CPU_HPP
