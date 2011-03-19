/*
 *  vli_number_gpu.hpp
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 19.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#include "vli_number_gpu.h"
#include "definition.h"
#include "GpuManager.h"



namespace vli
{
	/**
	 constructors 
	 */	
	template <class T>
	vli_gpu<T>::vli_gpu()
	{
		size_ = SIZE_BITS/(8*sizeof(T)); 
		gpu::check_error( cublasAlloc( size_, sizeof(T), (void**)&data_ ), __LINE__);	
	}		
	
	template <class T>
	vli_gpu<T>::vli_gpu(T num)
	{
		size_ = SIZE_BITS/(8*sizeof(T)); 
		gpu::check_error( cublasAlloc( size_, sizeof(T), (void**)&data_ ), __LINE__);	
		std::vector<T>  Array(size_,0);
		Array[0] = num;
		check_error(cublasSetVector(Array.size(), sizeof(T), &Array[0],1,p(), 1), __LINE__);
	}		

	template <class T>
	vli_gpu<T>::vli_gpu(vli_gpu<T> const& vli):size_(vli.size())
	{
		gpu::check_error(   cublasAlloc(size_,sizeof(T), (void**) &data_ ), __LINE__);
		cudaMemcpy( data_, vli.data_, size_*sizeof(T) , cudaMemcpyDeviceToDevice);			
	}
	
	template <class T>
	vli_gpu<T>::vli_gpu(vli_cpu<T> const& vli):size_(vli.size())
	{
		gpu::check_error(cublasAlloc(size_,sizeof(T), (void**) &data_ ), __LINE__);
		gpu::check_error(cublasSetVector(vli.size(), sizeof(T), &vli[0],1,p(), 1), __LINE__);
	}
	
	/**
	 destructors
	 */
	template <class T>
	vli_gpu<T>::~vli_gpu()
	{
		cublasFree(data_);
	}
	
	template <class T>
	vli_gpu<T> & vli_gpu<T>::operator = (vli_gpu<T> const &  vli)
	{
		swap(vli);
		return *this;
	}

	
	template <class T>
	void vli_gpu<T>::operator = (vli_cpu<T> const &  vli)
	{
		gpu::check_error(cublasSetVector(vli.size(), sizeof(T), &vli[0],1,p(), 1), __LINE__);
	}
	
	
	template <class T>
	void vli_gpu<T>::swap(vli_gpu& vli)
	{
		std::swap(data_, vli.data_);	
	}		
	
	template <class T>
	void vli_gpu<T>::copy_vli_to_cpu(vli::vli_cpu<T>& vli) 
	{
		check_error(cublasGetVector(vli.size(), sizeof(T), data_  ,1,&vli[0], 1), __LINE__);			
	}
	
	template <class T>
	vli_gpu<T>::operator vli::vli_cpu<T>()
	{
		vli::vli_cpu<T> vli(this->size_);
		copy_vli_to_cpu(vli);
		return vli;
	}
	
	template <class T>
	std::size_t vli_gpu<T>::size() const
	{
		return size_;
	}
	
	template <class T>
	inline const T* vli_gpu<T>::p () const
	{
		return data_;
	}
	
	template <class T>
	inline  T* vli_gpu<T>::p () 
	{
		return data_;
	}
	
	template <class T>
	std::ostream& operator<<(std::ostream& os, const vli_gpu<T> & vli_gpu)
	{
		int size = vli_gpu.size();		
		std::vector<T>  Array(size,0);
		
		cublasGetVector(vli_gpu.size(), sizeof(T), vli_gpu.p()  ,1,&Array[0], 1);
		
		for(typename std::vector<T>::const_iterator it = Array.begin();it !=Array.end();it++)
			os << *it<< " ";
		os<<std::endl;
		
		return os;		
		
	}
	
	
}