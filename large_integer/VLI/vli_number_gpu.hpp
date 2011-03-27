/*
 *  vli_number_gpu.hpp
 *  vli_cpu_gpu
 *
 *  Created by Tim Ewart on 19.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#include "definition.h"
#include "vli_number_gpu.h"
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
		gpu::check_error(cublasSetVector(Array.size(), sizeof(T), &Array[0],1,p(), 1), __LINE__);
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
	// should be possible to think
	template <class T>
	void vli_gpu<T>::operator = (vli_gpu<T> const &  vli)
	{
		std::cout << " srljkgs " << vli << std::endl;
		
		gpu::check_error( cudaMemcpy( this->p(), vli.p(), this->size()*sizeof(T) , cudaMemcpyDeviceToDevice), __LINE__);

	}
	/*
	template <class T>
	void vli_gpu<T>::operator = (vli_cpu<T> const &  vli)
	{
		gpu::check_error(cublasSetVector(vli.size(), sizeof(T), &vli[0],1,p(), 1), __LINE__);
	}
	*/
	template <class T>
	vli_gpu<T> & vli_gpu<T>::operator += (vli_gpu<T> const &  vli)
	{
		plus_assign(*this, vli);
		return *this;
	}
	
	template <class T>
	vli_gpu<T> & vli_gpu<T>::operator *= (vli_gpu<T> const &  vli)
	{
		multiply_assign(*this, vli);
		return *this;
	}
	
	template <class T>
	void vli_gpu<T>::swap(vli_gpu const& vli)
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
	T vli_gpu<T>::size() const
	{
		return size_;
	}
	
	template <class T>
	inline const T* vli_gpu<T>::p() const
	{
		return data_;
	}
	
	template <class T>
	inline  T* vli_gpu<T>::p() 
	{
		return data_;
	}
	
	
	template <class T>
	const vli_gpu<T> operator+ (vli_gpu<T> vli_gpu_a , vli_gpu<T> const & vli_gpu_b)
	{
		vli_gpu_a += vli_gpu_b;
	}
	
	template<class T>
	void plus_assign(vli_gpu<T>& vli_a, vli_gpu<T> const & vli_b)
	{
		addition_gpu(vli_a.p(),vli_b.p(),1,vli_b.size());
	}

	template <class T>
	const vli_gpu<T> operator* (vli_gpu<T> vli_gpu_a , vli_gpu<T> const & vli_gpu_b)
	{
		vli_gpu_a *= vli_gpu_b;
	}
	
	template<class T>
	void multiply_assign(vli_gpu<T>& vli_a, vli_gpu<T> const & vli_b)
	{
		multiply_gpu(vli_a.p(),vli_b.p(),1,vli_b.size());
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