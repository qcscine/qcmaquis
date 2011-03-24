/*
 *  vli_vector_cpu.hpp
 *  Untitled
 *
 *  Created by Tim Ewart on 23/03/11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#include "vli_vector_cpu.h"

namespace vli
{
	template<class T>
	vli_vector<T>::vli_vector(std::size_t size = 0):size_(size)
	{
		vector_vli_.resize(size);
		typename std::vector< vli_cpu<T> * >::iterator it;
		for(it = vector_vli_.begin() ; it != vector_vli_.end() ; it++)
		{
			*it = new vli_cpu<T>();
		}
	}
	
	template<class T>
	vli_vector<T>::~vli_vector()
	{
		vector_vli_.clear(); //Useless I think
	};


	template<class T>
	void vli_vector<T>::clean()
	{
		vector_vli_.clear(); //should be enough
	}
	
	template<class T>
	const std::size_t vli_vector<T>::GetSize() const
	{
		return size_;
	}
	
	/**
	 access operator
	 */
	template<class T>
	vli_cpu<T> const &  vli_vector<T>::operator[](std::size_t i ) const
	{
		return *(vector_vli_[i]);
	}
	
	template<class T>
	vli_cpu<T>  &  vli_vector<T>::operator[](std::size_t i )
	{
		return *(vector_vli_[i]);
	}
	
	template<class T>
	vli_vector<T>& vli_vector<T>::operator += (vli_vector<T> const &  vli)
	{
		plus_assign(*this, vli);
		return *this;
	}

	template<class T>
	vli_vector<T>& vli_vector<T>::operator *= (vli_vector<T> const &  vli)
	{
		multiply_assign(*this, vli);
		return *this;
	}	
	
	/**
	 multiply and addition operators, suite ...
	 */
	template <class T>
	const vli_vector<T> operator+(vli_vector<T> vli_a, vli_vector<T> const& vli_b)
	{
		vli_a += vli_b;
		return vli_a;
	}
																				   																				   																		
	template <class T>
	void plus_assign(vli_vector<T> & vli_a, vli_vector<T> const& vli_b )
	{
		/**
		  the vectors have the same size
		*/
		#pragma omp parallel for private(i)
		for(std::size_t i = 0; i < vli_b.GetSize(); i++)
		{
			addition_classic_cpu( vli_a[i], vli_b[i]);
		}
	}
	
	template <class T>
	const vli_vector<T> operator*(vli_vector<T> vli_a, vli_vector<T> const& vli_b)
	{
		vli_a *= vli_b;
		return vli_a;
	}
	
	template <class T>
	void multiply_assign(vli_vector<T> & vli_a, vli_vector<T> const& vli_b )
	{
		/**
		 the vectors have the same size
		 */
		#pragma omp parallel for private(i)
		for(std::size_t i = 0; i < vli_b.GetSize(); i++)
		{
			multiplication_classic_cpu( vli_a[i], vli_b[i]);
		}
	}
	
	
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