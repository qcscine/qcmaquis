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
		cublasInit();
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
		cublasShutdown();
	};
	
	
	cudaDeviceProp gpu_manager::GetDeviceProperties() const
	{
		return deviceProp_;
	}
	
	
	
	void check_error(cublasStatus const& stat, std::size_t line)
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
	
}

