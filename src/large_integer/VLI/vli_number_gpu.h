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

#include <boost/lexical_cast.hpp>

namespace vli
{
	
	
	template<class T>
	class vli_gpu : public bvli
	{
	public:
		vli_gpu()
		{
			size_ = SIZE_BITS/(8*sizeof(T)); // size of the vli, understand number of box
			check_error( cublasAlloc( size_, sizeof(T), (void**)&data_ ), __LINE__);	
		}		
				
		vli_gpu(T num)
		{
			size_ = SIZE_BITS/(8*sizeof(T)); // size of the vli, understand number of box
			check_error( cublasAlloc( size_, sizeof(T), (void**)&data_ ), __LINE__);	
			std::vector<T>  Array(size_,0);
			Array[0] = num;
			check_error(cublasSetVector(Array.size(), sizeof(T), &Array[0],1,p(), 1), __LINE__);
		}		
		
		~vli_gpu()
		{
			cublasFree(data_);
		}
		
		std::size_t size() const
		{
			return size_;
		}
		
		inline const T* p () const
		{
			return data_;
		}
		
		inline  T* p () 
		{
			return data_;
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
					//std::cout << "CUBLAS_STATUS_SUCCESS" + error << std::endl;
					break;
			}
			
		}
	
	private:
		T* data_;
		int size_;
		
	};
	
	
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
#endif
