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

#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/matrix_interface.hpp"
#include "p_dense_matrix/resizable_matrix_interface.hpp"
#include "p_dense_matrix/dense_matrix_algorithms.h"
#include "p_dense_matrix/matrix_algorithms.hpp"

#include <boost/timer.hpp>

#include <omp.h>
#include "p_dense_matrix/dense_matrix_blas.hpp"


int main(int   argc, char * argv[])
{
	//gpu::Simu Simulation;
	
	srand(0);
	int NUM = atoi(argv[1]);
	
	blas::p_dense_matrix<float> A(2*NUM,3*NUM);
	blas::p_dense_matrix<float> B(3*NUM,2*NUM);
	blas::p_dense_matrix<float> C(2*NUM,3*NUM,0);
	
	for(int i=0; i< 2*NUM ;i++  )
	{
		for(int j=0; j< 3*NUM ; j++ )
		{
			A(i,j) = rand();
		}
	}

	cout << A << endl;
	
	for(int i=0; i< 3*NUM ;i++  )
	{
		for(int j=0; j<2*NUM ; j++ )
		{
				B(i,j) = rand();
		}
	}

	cout << B << endl;
	
	
	struct timeval tp;
	gettimeofday( &tp, NULL );
		
	double sec      = static_cast<double>( tp.tv_sec );
	double usec = static_cast<double>( tp.tv_usec )/1E6;
	double start = sec + usec;
	
	C = matrix_matrix_multiply(A,B);
	
	cout << C << endl;
	
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
	
	out.open("resultatcpu.txt",std::ios::app);
	out << C <<endl;
	out.close();
}
