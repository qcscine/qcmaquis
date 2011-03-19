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
		
		vli_gpu<T> & operator = (vli_gpu<T> const &  vli);
		void operator = (vli_cpu<T> const &  vli);

		void swap(vli_gpu& vli);
		void copy_vli_to_cpu(vli::vli_cpu<T>& vli); 
		operator vli::vli_cpu<T>();
		
		
		std::size_t size() const ;
		inline const T* p () const ;
		inline  T* p (); 

	private:
		T* data_;
		int size_;
		
	};
	
	template <class T>
	std::ostream& operator<<(std::ostream& os, const vli_gpu<T> & vli_gpu);
	

	
}
#endif
