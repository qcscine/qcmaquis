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

#include <sys/time.h> 
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;

#define NDEBUG

#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"

#include <boost/timer.hpp>

#include <omp.h>


#ifdef __bGPU__
#include "dense_matrix/gpu/matrix_gpu.h"
#include "dense_matrix/gpu/dense_matrix_gpu.h"
#include "dense_matrix/gpu/matrix_gpu_functions.hpp"
#else
#include "dense_matrix/dense_matrix_blas.hpp"
#endif


int main(int   argc, char * argv[])
{
//	gpu::Simu Simulation;
	
	srand(0);
	int NUM = atoi(argv[1]);
	
	blas::dense_matrix<float> A(NUM,NUM);
	blas::dense_matrix<float> B(NUM,NUM);
	blas::dense_matrix<float> C(NUM,NUM,0);
	
	for(int i=0; i++ ; i< NUM)
	{
		for(int j=0; j++ ; j< NUM)
		{
			A(i,j) = rand();
			B(i,j) = rand();
		}
	}
	struct timeval tp;
	gettimeofday( &tp, NULL );
	double sec      = static_cast<double>( tp.tv_sec );
	double usec = static_cast<double>( tp.tv_usec )/1E6;
	double start = sec + usec;
	
	C = matrix_matrix_multiply(A,B);
	
	
	gettimeofday( &tp, NULL );
	sec = static_cast<double>( tp.tv_sec );
	usec = static_cast<double>( tp.tv_usec )/1E6;
	double end = sec + usec;
	
	
	double time = end - start;
	
	std::ofstream out;
	out.open("time.txt",std::ios::app);
	out << std::setprecision (6) << time << " " <<  NUM << " " << omp_get_max_threads() <<endl;
	out.close();
	
	cout << std::setprecision (6) << time << " " <<  NUM << " " << omp_get_max_threads() <<endl;
	
}
