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
#include <vector>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include "vector.hpp"


namespace gpu
{


template<class T>
class vector_gpu
{
public:
    typedef T   value_type;
    typedef std::size_t  size_type;


	vector_gpu(size_type size):size_(size) 
	{
		check_error(cublasAlloc( size, sizeof(T), (void**)&p_ ), __LINE__);
	}
	
	vector_gpu(size_type size, T value):size_(size)
	{
        check_error( cublasAlloc( size, sizeof(T), (void**)&p_ ), __LINE__);
	
		std::vector<T>  Array(size, value);		
		check_error( cublasSetVector(size,sizeof(T),&Array[0],1,p(),1), __LINE__);
	}

    vector_gpu(vector_gpu const& v)
        : size_(v.size_)
    {
        check_error( cublasAlloc( size_, sizeof(T), (void**)&p_), __LINE__);
		check_error( cudaMemcpy( p_, v.p_, size_*sizeof(T) , cudaMemcpyDeviceToDevice), __LINE__);
    }
		
	template<class MemoryBlock>
	vector_gpu(maquis::types::vector<T, MemoryBlock> const & Vector_cpu):size_(Vector_cpu.size())
	{
		check_error( cublasAlloc( size_, sizeof(T), (void**)&p_ ), __LINE__);	
		check_error( cublasSetVector(size_,sizeof(T), &Vector_cpu(0),1,p_,1), __LINE__);
	};

	template <typename MemoryBlock>
	void copy_vector_to_cpu(maquis::types::vector<T,MemoryBlock>& v_cpu) const
	{
		check_error( cublasGetVector(size_, sizeof(T),p(),1,&v_cpu(0),1), __LINE__);			
	}
	
	template <typename MemoryBlock>
	operator maquis::types::vector<T,MemoryBlock>()
	{
		maquis::types::vector<T,MemoryBlock> v(size_);
		copy_vector_to_cpu(v);
		return v;
	}
	
	~vector_gpu()
	{
		cublasFree(p_);
	}

    void swap(vector_gpu& v)
    {
        std::swap(size, v.size_);
        std::swap(p_, v.p_);
    }

    friend void swap(vector_gpu& v1, vector_gpu& v2)
    {
        v1.swap(v2);
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

    vector_gpu& operator = (vector_gpu v)
    {
        swap(v);
        return *this;
    }

	inline void check_error(cublasStatus const& stat, unsigned int line)
	{
		switch (stat) 
		{
			case CUBLAS_STATUS_NOT_INITIALIZED:
                throw(std::runtime_error("CUBLAS_STATUS_NOT_INITIALIZED in " + boost::lexical_cast<std::string>(__FILE__) + boost::lexical_cast<std::string>(line) ));
				break;
				
			case CUBLAS_STATUS_MAPPING_ERROR:
                throw(std::runtime_error("CUBLAS_STATUS_MAPPING_ERROR in " + boost::lexical_cast<std::string>(__FILE__) + boost::lexical_cast<std::string>(line) ));
				break;
				
			case CUBLAS_STATUS_INVALID_VALUE:
                throw(std::runtime_error("CUBLAS_STATUS_INVALID_VALUE in " + boost::lexical_cast<std::string>(__FILE__) + boost::lexical_cast<std::string>(line) ));
				break;	
				
			default:
				//maquis::cout << "CUBLAS_STATUS_SUCCESS" + error << std::endl;
                break;
		}
		
	}
	
private:
	T* p_;
	size_type size_;
};

}
#endif
