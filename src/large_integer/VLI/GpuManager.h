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
	
	/**
	 checking error function
	 */
	void check_error(cublasStatus const& stat, std::size_t line);

}

#endif
