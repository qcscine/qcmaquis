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
#include <cublas.h>

#include "boost/lexical_cast.hpp"
#include <shrUtils.h>
#include <cuda_runtime_api.h>
#include <stdexcept>

namespace gpu
{
	class gpu_manager
	{
	private:
		gpu_manager(int device = 0):device_(device)
		{
			cuInit(device_);
			cublasInit();
			cudaGetDeviceProperties(&deviceProp_, device_);
		};			
		gpu_manager(gpu_manager const&);
		gpu_manager& operator=(gpu_manager const&);
	public:
		~gpu_manager()
		{
			destructor();
		}
		
		static gpu_manager& instance();
		
		cudaDeviceProp const GetDeviceProperties()
		{
			return deviceProp_;
		}
		
		void constructor();
		void destructor();
		
		
			
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
	

}

#endif