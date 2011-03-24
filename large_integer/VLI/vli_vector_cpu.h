/*
 *  vli_vector_cpu.h
 *  Untitled
 *
 *  Created by Tim Ewart on 23/03/11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef __VLI_VECTOR_
#define __VLI_VECTOR_

#include "definition.h"
#include "vli_number_cpu.h"

namespace vli
{
template <class T>
class vli_vector
{
public:

	vli_vector(std::size_t size);
	~vli_vector();
	void clean();
	std::size_t const GetSize() const;

	/**
	 access operator
	*/
	vli_cpu<T>const& operator[](std::size_t i) const;
	vli_cpu<T>& operator[](std::size_t i) ;	
	/**
	 multiply and addition operators
	 */															
	vli_vector<T>& operator += (vli_vector<T> const& vli); 
	vli_vector<T>& operator *= (vli_vector<T> const& vli); 
	
private:
	std::size_t size_;
	std::vector<vli_cpu<T>* > vector_vli_;
};
	
	/**
	 multiply and addition operators, suite ...
	 */
	template <class T>
	const vli_vector<T> operator+(vli_vector<T> vli_a, vli_vector<T> const& vli_b);
	
	template <class T>
	void plus_assign(vli_vector<T> & vli_a, vli_vector<T> const& vli_b );
	
	template <class T>
	const vli_vector<T> operator*(vli_vector<T> vli_a, vli_vector<T> const& vli_b);
	
	template <class T>
	void multiply_assign(vli_vector<T> & vli_a, vli_vector<T> const& vli_b );
	
	/**
	 stream 
	 */
	template<typename T>
	std::ostream& operator<< (std::ostream& os,  vli_vector<T> & vli);
	
	

}

#endif