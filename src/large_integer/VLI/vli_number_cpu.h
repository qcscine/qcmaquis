/*
 *  vli.h
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */


#ifndef __VLI__
#define __VLI__

#include "definition.h"

namespace vli
{
	
	/**
	template forward declaration 
	*/
	template<class T>
	class vli_gpu;

	template<class T>
	class vli_cpu 
	{
	public:
		
		/**
		constructors 
		*/
		vli_cpu();
		vli_cpu(T num);
		
		/**
		destructor 
		*/
		~vli_cpu();
		
		/**
		 copy constructors, CPU to CPU and GPU to CPU
		*/
		vli_cpu<T> & operator= (vli_cpu<T>  vli);
		void swap(vli_cpu& vli);
		operator vli::vli_gpu<T>();
		void copy_vli_to_gpu(vli::vli_gpu<T>& vli) const;
		
		/**
		logistics operators 
		*/
		size_int size() const; 
		T &operator[](T i); 
		T const & operator[](T i) const;
		typename std::vector<T>::const_iterator begin();
		typename std::vector<T>::const_iterator end();

		/**
		multiply and addition operators
		*/
		vli_cpu<T>& operator += (vli_cpu<T> const& vli); 
		vli_cpu<T>& operator *= (vli_cpu<T> const& vli); 
		
	private:
		std::vector<T> data_;
		size_int size_;
		
	};

	/**
	 multiply and addition operators, suite ...
	 */
	template <class T>
	const vli_cpu<T> operator+(vli_cpu<T> vli_a, vli_cpu<T> const& vli_b);
	
	template <class T>
	void plus_assign(vli_cpu<T> & vli_a, vli_cpu<T> const& vli_b );
	
	template <class T>
	const vli_cpu<T> operator*(vli_cpu<T> vli_a, vli_cpu<T> const& vli_b);
	
	template <class T>
	void multiply_assign(vli_cpu<T> & vli_a, vli_cpu<T> const& vli_b );
	
	/**
	stream 
	*/
	template<typename T>
	std::ostream& operator<< (std::ostream& os,  vli_cpu<T> & vli);
	 
}

#endif
