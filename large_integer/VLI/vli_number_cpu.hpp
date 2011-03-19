/*
 *  vli_num_cpu.hpp
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 19.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#include "vli_number_cpu.h"
#include "kernels_cpu.h"
#include "GpuManager.h"


/**
 TO DO !
 Introduce a kind of hook, to make the choise between the classic, Avizienis, gmp or something else insign plus_assign, etc ..
 Andreas proposes something ^^'
 */

namespace vli
{
	/**
	 constructors 
	 */
	template <class T>	
	vli_cpu<T>::vli_cpu()
	{
		size_ = SIZE_BITS/(8*sizeof(T)); 
		data_.resize(size_);
	}	
	
	template <class T>		
	vli_cpu<T>::vli_cpu(T num)
	{
		size_ = SIZE_BITS/(8*sizeof(T)); 
		data_.resize(size_);
		data_[0] = num;
	}		
	
	/**
	 destructor 
	 */
	template <class T>
	vli_cpu<T>::~vli_cpu()
	{
		data_.erase(data_.begin(),data_.end());
	}
	
	/**
	 copy constructors, CPU to CPU and GPU to CPU
	 */
	template <class T>
	vli_cpu<T>& vli_cpu<T>::operator=(vli_cpu<T> vli)
	{
		swap(vli);
		return *this;
	}
	
	template <class T>
	void vli_cpu<T>::swap(vli_cpu& vli)
	{
		std::swap(data_, vli.data_);	
	}
	
	template <class T>
	vli_cpu<T>::operator vli::vli_gpu<T>()
	{
		vli_cpu<T> vli(this->size_);
		copy_vli_to_gpu(vli);
		return vli;
	}
	
	template <class T>
	void vli_cpu<T>::copy_vli_to_gpu(vli::vli_gpu<T>& vli) const
	{
		gpu::check_error(cublasSetVector(this->size(), sizeof(T), &this->data_[0]  ,1,vli.p(), 1), __LINE__);			
	}
	
	/**
	Logistics 
	*/
	template <class T>
	size_int vli_cpu<T>::size() const
	{
		return size_;
	}

	template <class T>
	T& vli_cpu<T>::operator[](T i)
	{
		return data_[i];
	}
	
	template <class T>
	T const & vli_cpu<T>::operator[](T i) const
	{
		return data_[i];
	}

	template <class T>
	typename std::vector<T>::const_iterator vli_cpu<T>::begin()
	{
		typename std::vector<T>::const_iterator it= data_.begin();
		return it;
	}
	
	template <class T>
	typename std::vector<T>::const_iterator vli_cpu<T>::end()
	{
		typename std::vector<T>::const_iterator it= data_.end();
		return it;
	}
	

	/**
	additions and multiplications 
	*/
	template <class T>
    const vli_cpu<T> operator+ (vli_cpu<T> a, vli_cpu<T> const& b)
    {
        a += b;
        return a;
    }
	
	template <class T>
	vli_cpu<T>&  vli_cpu<T>::operator += (vli_cpu<T> const& vli)
	{
		plus_assign(*this,vli);
		return *this;
	}

	template <class T>
	void plus_assign(vli_cpu<T> & vli_a, vli_cpu<T> const& vli_b)
	{
		addition_classic_cpu(vli_a,vli_b);
	}
	
	template <class T>
	const vli_cpu<T> operator* (vli_cpu<T> a, vli_cpu<T> const& b)
	{
		a *= b;
		return a;
	}
	
	template <class T>
	vli_cpu<T>& vli_cpu<T>::operator *= (vli_cpu<T> const& vli)
	{
		multiply_assign(*this , vli);
		return *this;
	}
	
	template <class T>
	void multiply_assign(vli_cpu<T> & vli_a, vli_cpu<T> const& vli_b)
	{
		multiplication_classic_cpu(vli_a, vli_b);
	}
	
	/**
	stream 
	*/
	template<typename T>
	std::ostream& operator<< (std::ostream& os,  vli_cpu<T> & vli)
	{
		for(typename std::vector<T>::const_iterator it=vli.begin(); it != vli.end(); it++)
			std::cout << *it << " " ;
			os<<std::endl; 
			
			return os;
	}
		
	
}