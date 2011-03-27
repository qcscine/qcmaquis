/*
 *  vli.h
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef __VLI_GPU__
#define __VLI_GPU__


#include <iostream>
#include <stdexcept>
#include <assert.h>

#include "boost/lexical_cast.hpp"


#include "definition.h"
#include "vli_number_cpu.h"
#include "GpuManager.h"

namespace vli
{

		
	template<class T>
	class vli_gpu 
	{
	public:
		/**
		constructors 
		*/
		vli_gpu();
		vli_gpu(T num);
		vli_gpu(vli_gpu<T> const& vli);
		vli_gpu(vli_cpu<T> const& vli);

		
		/**
		destructors 
		*/
		~vli_gpu();
		//think on the copy swap, I think it is totaly useless for gpu
	//	vli_gpu<T> & operator = (vli_gpu<T> const &  vli);
		/**
		due to gpu architecture the overload is weird
		*/
	//	void operator = (vli_cpu<T> const &  vli);
		void operator = (vli_gpu<T> const &  vli);

		void swap(vli_gpu const& vli);
		void copy_vli_to_cpu(vli::vli_cpu<T>& vli); 
		operator vli::vli_cpu<T>();
		
		T size() const ;
		inline  const T* p() const ;
		inline  T* p(); 

		/**
		 multiply and addition operators
		 */
		vli_gpu<T>& operator += (vli_gpu<T> const& vli); 
		vli_gpu<T>& operator *= (vli_gpu<T> const& vli); 
		
	private:
		T* data_;
		T size_;
		
	};
	
	/**
	 multiply and addition operators, suite ...
	 */
	template <class T>
	const vli_gpu<T> operator+(vli_gpu<T> vli_a, vli_gpu<T> const& vli_b);
	
	template <class T>
	void plus_assign(vli_gpu<T> & vli_a, vli_gpu<T> const& vli_b );
	
	template <class T>
	const vli_gpu<T> operator*(vli_gpu<T> vli_a, vli_gpu<T> const& vli_b);
	
	template <class T>
	void multiply_assign(vli_gpu<T> & vli_a, vli_gpu<T> const& vli_b );
	
	template <class T>
	std::ostream& operator<<(std::ostream& os, const vli_gpu<T> & vli_gpu);
	
}
#endif
