/*
 *  matrix_gpu.cpp
 *  XCODE_MAGMA
 *
 *  Created by Tim Ewart on 29.11.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */

#ifndef __MATRIX_GPU_FUNCTIONS__
#define __MATRIX_GPU_FUNCTIONS__

#include <iostream>
#include <cassert>

#include <vector>


#include "cublas.h"

#include "matrix_gpu.h"
#include "timings.h"

#include "definition.h"


/*
My GT 330 does not support double so I develop, and debug with float.
run on CSCS with double, moreover we must respect the f77 philosophy.
*/


namespace gpu 
{
	

	
template<class T>
void addition_classic_gpu(vli::vli_matrix<T> const & lhs,vli::vli_matrix<T> const & rhs,vli::vli_matrix<T> & result_cpu)
{
	matrix_gpu<T> lhs_GPU(lhs);
	matrix_gpu<T> rhs_GPU(rhs);
	matrix_gpu<T> result_cpu_GPU(result_cpu);
	
	int num_integer  = lhs_GPU.num_columns();
	int ld           = lhs_GPU.num_rows();
		
	addition(lhs_GPU.p(),rhs_GPU.p(),result_cpu_GPU.p(), num_integer, ld);
		
	std::cout << result_cpu_GPU << std::endl;

}


template <class T>
std::ostream& operator<< (std::ostream& os, const  matrix_gpu<T> & Matrix_gpu)
{
	size_type size1 = Matrix_gpu.size1();
	size_type size2 = Matrix_gpu.size2();			
	size_type reserved_size1 = Matrix_gpu.size1();
	
	std::vector<T>  Array(size1*size2);
		
	cublasGetMatrix(size1,size2,sizeof(T),Matrix_gpu.p(),size1,&Array[0],size1);
		
	for(size_type i=0; i< size1; ++i)
	{
		for(size_type j=0; j < size2; ++j)
			os << Array[i+j*reserved_size1] << " ";
		os << std::endl;
	}
	
	return os;
}

	
	
	
}

#endif

