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

#include <shrUtils.h>
#include <cuda_runtime_api.h>

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
	

}

#endif