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

#define NDEBUG

#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"

#include <omp.h>

#include "dense_matrix/dense_matrix_blas.hpp"
#include "dense_matrix/gpu/matrix_gpu.h"
#include "dense_matrix/gpu/matrix_gpu_functions.hpp"



template <class T>
inline std::string to_string (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}


int main( int argc, char* argv[])
{
	
	gpu::Simu Simulation;
	srand(0);
	
	double verif = 0 ;
	double verif2 = 0 ;

	
	for (int NUM = 3 ; NUM <= 15000 ; NUM *= 2)
	{
		
		blas::dense_matrix<double> A(NUM,NUM);
		blas::dense_matrix<double> B(NUM,NUM);
		blas::dense_matrix<double> C_CPU(NUM,NUM,0);
		blas::dense_matrix<double> C_GPU(NUM,NUM,0);
		
		for(int i=0; i< NUM ;i++  )
		{
			for(int j=0; j< NUM ; j++ )
			{
				A(i,j) = rand();
			}
		}
		
		for(int i=0; i< NUM ;i++  )
		{
			for(int j=0; j< NUM ; j++ )
			{
				B(i,j) = rand();
			}
		}
		
		string size_n = to_string(NUM);
		string name_CPU = "DGEMM CPU ";
		string name_CPU_GPU = "DGEMM GPU ";
		
		Timer temps_CPU(name_CPU += size_n) ;
		
		temps_CPU.begin();	
		C_CPU =  matrix_matrix_multiply(A,B);
		temps_CPU.end();
		
		double N = static_cast<double>(NUM);
		double GFLOP_CPU = N*N*N/temps_CPU.GetTime();
		
		
		
		Timer temps_CPU_GPU(name_CPU_GPU += size_n) ;
		
		temps_CPU_GPU.begin();	
		gpu::matrix_matrix_multiply(A,B,C_GPU);
		temps_CPU_GPU.end();	
	
		cout <<temps_CPU_GPU.GetTime() << endl;	

	
		double GFLOP_GPU =2*N*N*N/temps_CPU_GPU.GetTime();
		
		double split = GFLOP_GPU/(GFLOP_GPU+GFLOP_CPU);
		double split2 =  N*sqrt(GFLOP_GPU/(GFLOP_GPU+GFLOP_CPU)); 
		
		cout << " GFLOP_GPU   "  << GFLOP_GPU << " GFLOP_CPU " <<  GFLOP_CPU << endl; 
		cout << " split  " << split << endl;
	/*	
		for(int i=0; i< NUM ;i++  )
		{
			for(int j=0; j< NUM ; j++ )
			{
				verif = (C_GPU(i,j) -  C_CPU(i,j) ) * (C_GPU(i,j) -  C_CPU(i,j) ) / C_CPU(i,j);
			}
		}
	*/	
		
		std::ofstream out;
		out.open("time.txt",std::ios::app);
		out << std::setprecision (6) <<  NUM << " " << temps_CPU.GetTime() << " " << temps_CPU_GPU.GetTime() << "   " << split << " " << split2 << " "  << omp_get_max_threads() <<endl;
		out.close();
		
		
	}	
	return 0;
}

