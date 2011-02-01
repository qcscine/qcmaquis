/*
 *  transposecpu_gpu.cpp
 *  dmrg
 *
 *  Created by Tim Ewart on 28.01.11.
 *  Copyright 2011 Université de Genève. All rights reserved.
 *
 */

/*
 *  dgemmcpugpu.cpp
 *  dmrg
 *
 *  Created by Tim Ewart on 25/12/10.
 *  Copyright 2010 Université de Genève. All rights reserved.
 *
 */
#include <cmath>
#include <iterator>
#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>


#include <sys/time.h> 
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::string;



#include "utils/timings.h"


#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"

#include <omp.h>

#include "dense_matrix/dense_matrix_blas.hpp"
#include "dense_matrix/gpu/matrix_gpu.h"
#include "dense_matrix/gpu/matrix_gpu_functions.hpp"
#include "dense_matrix/gpu/GpuManager.h"

#include "dense_matrix/gpu/list_functions.h"

using namespace gpu;



#define NUM 17

#define NUM2 17



int main( int argc, char* argv[])
{
	gpu::gpu_manager*  GPU;
	GPU->instance();
	
	blas::dense_matrix<float> A(NUM,NUM2);
	blas::dense_matrix<float> B(NUM2,NUM);

	
	for(int i=0; i < NUM ; i ++)
	{
		for(int j=0; j < NUM2 ; j ++)
		{	
			A(i,j) = j;
		}
	}
		
		
//	std::cout << A << std::endl;
	
	gpu::matrix_gpu<float> A_GPU(A);
	gpu::matrix_gpu<float> B_GPU(B);
	
//	transpose(B_GPU.p(),A_GPU.p(),A_GPU.num_rows(),A_GPU.num_columns(),A_GPU.ld());
	std::cout << " out avant swap : " << std::endl;
	std::cout << A_GPU << std::endl;
	
//	swap_rows(A_GPU.p(),A_GPU.num_rows(),A_GPU.num_columns(),A_GPU.ld(),3,80);
	swap_columns(A_GPU.p(),A_GPU.num_rows(),B_GPU.num_columns(),B_GPU.ld(),3,8);
	
//	std::cout << " in : " << std::endl;
//	std::cout << A_GPU << std::endl;

	std::cout << " out apres swap : " << std::endl;
	std::cout << A_GPU << std::endl;
	
	
	GPU->instance().destructor();	
};


