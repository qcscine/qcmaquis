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

namespace gpu
{
	class gpu_manager
	{
	private:
		gpu_manager()
		{
			cuInit(0);
			cublasInit();
		};			
		gpu_manager(gpu_manager const&);
		gpu_manager& operator=(gpu_manager const&);
	public:
		static gpu_manager& instance();
		void constructor();
		void destructor();
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