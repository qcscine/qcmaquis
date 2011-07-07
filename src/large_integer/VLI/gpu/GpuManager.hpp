/*
 *  GpuManager.cpp
 *  Untitled
 *
 *  Created by Tim Ewart on 24/03/11.
 *  Copyright 2011 Univerisité de Geneve . All rights reserved.
 *
 */

#include "GpuManager.h"

namespace gpu 
{
	gpu_manager::gpu_manager(int device = 0):device_(device)
	{
		cuInit(device_);
		cudaGetDeviceProperties(&deviceProp_, device_);
	};			

	gpu_manager::~gpu_manager()
	{
		destructor();
	}
	gpu_manager& gpu_manager::instance()
	{
		static gpu_manager* singleton = NULL;
		if (!singleton)
		{
			singleton = new gpu_manager();
		}
		return *singleton;
	};
	
	void gpu_manager::destructor()
	{
	};
	
	
	cudaDeviceProp gpu_manager::GetDeviceProperties() const
	{
		return deviceProp_;
	}
	
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
				//std::cout << "CUBLAS_STATUS_SUCCESS" + error << std::endl;
				break;
		}
	}
	
}

