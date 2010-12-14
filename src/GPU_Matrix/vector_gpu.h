/*
 *  vector_gpu.h
 *
 *  Created by Tim Ewart and Alex Kosemkov on 30.11.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */

#ifndef __VECTOR_GPU__
#define __VECTOR_GPU__

#include "cassert"

namespace gpu
{

typedef std::size_t  size_type;

template<class T>
class vector_gpu
{
public:
	vector_gpu(size_type size):size_(size) 
	{
		cublasAlloc( size, sizeof(T), (void**)&p_ );
	}
	
	vector_gpu(size_type size, T value):size_(size)
	{
		cublasAlloc( size, sizeof(T), (void**)&p_ );
		assert(true == CheckError(" cudaMalloc constructor vector"));
	
		std::vector<T>  Array(size, value);		
		cublasSetVector(size,sizeof(T),&Array[0],1,p(),1);
		assert(true == CheckError(" cublasSetVector constructor matrix"));
	}
		
	~vector_gpu()
	{
		cublasFree(p_);
	}
	
	
	T& operator()(const size_type i)
	{
		assert((i < size_));
		return p_[i];
	}

	inline const T* p () const
	{
		return p_;
	}
	
	inline  T* p () 
	{
		return p_;
	}
		
	//STL compatibility 
	inline const size_type size() const
	{
		return size_;
	}
	
	bool CheckError(std::string error)
	{
		switch (stat_) 
		{
			case CUBLAS_STATUS_NOT_INITIALIZED:
				std::cout << "CUBLAS_STATUS_NOT_INITIALIZED" + error << std::endl;
				return false;
				break;
				
			case CUBLAS_STATUS_MAPPING_ERROR:
				std::cout << "CUBLAS_STATUS_MAPPING_ERROR" + error << std::endl;
				return false;
				break;
				
			case CUBLAS_STATUS_INVALID_VALUE:
				std::cout << "CUBLAS_STATUS_INVALID_VALUE" + error << std::endl;
				return false;
				break;	
				
			default:
				std::cout << "CUBLAS_STATUS_SUCCESS" + error << std::endl;
				return true;
				break;
		}
		
	}
	
private:
	T* p_;
	size_type size_;
	cublasStatus stat_;
	
};

}
#endif