/*
 *  GpuManager.h
 *  dmrg
 *
 *  Created by Tim Ewart on 27.01.11.
 *  Copyright 2011 Université de Genève. All rights reserved.
 *
 */


#ifndef __GPU_MANAGER__
#define __GPU_MANAGER__

// utilities and system includes
#include <cuda.h>

#include "boost/lexical_cast.hpp"
#include <cuda_runtime_api.h>
#include <stdexcept>

namespace gpu
{
	class gpu_manager
	{
	private:
		gpu_manager(int device);
		gpu_manager(gpu_manager const&);
		gpu_manager& operator=(gpu_manager const&);
	public:
		~gpu_manager();
		static gpu_manager& instance();
		cudaDeviceProp GetDeviceProperties() const;		
		void constructor();
		void destructor();
        unsigned int GetmaxThreadsPerBlock();        
	private: 
		
		/**
			All fields of the devices are in
		*/
		cudaDeviceProp deviceProp_;	
		/**
			num of the device, will be link to ambient to avoid two init on the device 0
		*/
		int device_;
	};
	void cu_check_error(cudaError_t  const& err, std::size_t line)
	{
		switch (err) 
		{
			case cudaErrorInvalidDevicePointer:
				throw(std::runtime_error("cudaErrorInvalidDevicePointer, in " + boost::lexical_cast<std::string>(__FILE__) + boost::lexical_cast<std::string>(line) ));
				break;
				
			case cudaErrorMemoryAllocation:
				throw(std::runtime_error("cudaErrorMemoryAllocation in " + boost::lexical_cast<std::string>(__FILE__) + boost::lexical_cast<std::string>(line) ));
				break;
				
			case cudaErrorInitializationError:
				throw(std::runtime_error("cudaErrorInitializationError in " + boost::lexical_cast<std::string>(__FILE__) + boost::lexical_cast<std::string>(line) ));
				break;	
				
			case cudaErrorInvalidMemcpyDirection:	
				throw(std::runtime_error("cudaErrorInvalidMemcpyDirection in " + boost::lexical_cast<std::string>(__FILE__) + boost::lexical_cast<std::string>(line) ));
				break;	
			
			case  cudaErrorInvalidValue:
				throw(std::runtime_error("cudaErrorInvalidValue in " + boost::lexical_cast<std::string>(__FILE__) + boost::lexical_cast<std::string>(line) ));
				break;	
				
				
			default:
				//std::cout << "CUBLAS_STATUS_SUCCESS" << std::endl;
				break;
		}
	}
    
    void cu_check_error(cudaError_t  const& err, std::size_t line);
    
}


#endif
