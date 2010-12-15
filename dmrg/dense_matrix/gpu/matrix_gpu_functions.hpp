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
#include "cula.h"

#include "dense_matrix/gpu/matrix_gpu.h"
#include "dense_matrix/gpu/vector_gpu.h"

/*
My GT 330 does not support double so I develop, and debug with float.
run on CSCS with double, moreover we must respect the f77 philosophy.
*/

/*----------------------------- multiply -----------------------------------------*/

namespace gpu 
{

template<>
void matrix_gpu<float>::multiplies_assign(matrix_gpu<float> const& Matrix_right)
{
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = this->num_rows();
	size_type num_columns = Matrix_right.num_columns();
	size_type K =  this->num_columns();
	
	cublasSgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_columns, K, 1, this->p(), this->ld(), Matrix_right.p(), Matrix_right.ld(), 0, this->p_, this->ld());	
	
}
	
template<>
void gemm(matrix_gpu<float> const& A, matrix_gpu<float> const& B, matrix_gpu<float> & C)
{
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = A.num_rows();
	size_type num_columns = B.num_columns();
	size_type K =  A.num_columns();
	cublasStatus s; 
	cublasSgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_columns, K, 1, A.p(), A.ld(), B.p(), B.ld(), 0, C.p(), C.ld());	
	if (s != CUBLAS_STATUS_SUCCESS )
	{
		std::cout << " Sgemm " <<cublasGetError() << std::endl;
	}

}
	
template<>
void matrix_gpu<double>::multiplies_assign(matrix_gpu<double> const& Matrix_right)
{
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = this->num_rows();
	size_type num_columns = Matrix_right.num_columns();
	size_type K =  this->num_columns();
	
	cublasDgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_columns, K, 1, this->p(), this->ld(), Matrix_right.p(), Matrix_right.ld(), 0, this->p_, this->ld());	
	
}

template<>
matrix_gpu<float> operator* ( matrix_gpu<float> const & Matrix_left,  matrix_gpu<float>const & Matrix_right)
{
	assert(Matrix_left.num_rows() == Matrix_right.num_columns());
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_columns = Matrix_right.num_columns();
	size_type K =  Matrix_left.num_columns();
	
	matrix_gpu<float> Result(num_rows,num_columns);
	
	cublasSgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_columns, K, 1, Matrix_left.p(), Matrix_left.ld(), Matrix_right.p(), Matrix_right.ld(), 0, Result.p(), Result.ld());
	
	return Result;
};
	
template<>
matrix_gpu<float> matrix_matrix_multiply(matrix_gpu<float> const & Matrix_left,matrix_gpu<float> const & Matrix_right)
{
//	assert(Matrix_left.num_rows() == Matrix_right.num_columns());
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_columns = Matrix_right.num_columns();
	size_type K =  Matrix_left.num_columns();
	
	matrix_gpu<float> Result(num_rows,num_columns);
	
	cublasSgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_columns, K, 1, Matrix_left.p(), Matrix_left.ld(), Matrix_right.p(), Matrix_right.ld(), 0, Result.p(), Result.ld());
	
	return Result;
}



template<>
matrix_gpu<double> operator * ( matrix_gpu<double> const & Matrix_left,  matrix_gpu<double> const & Matrix_right)
{
	assert(Matrix_left.num_rows() == Matrix_right.num_columns());
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_columns = Matrix_right.num_columns();
	size_type K =  Matrix_left.num_columns();
	
	matrix_gpu<double> Result(num_rows,num_columns);
	
	cublasDgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_columns, K, 1, Matrix_left.p(), Matrix_left.ld(), Matrix_right.p(), Matrix_right.ld(), 0, Result.p(), Result.ld());
	
	return Result;
};

/*----------------------------- addition -----------------------------------------*/

template<>
void matrix_gpu<float>::plus_assign(matrix_gpu<float> const& Matrix_right)
{
	assert( (this->num_rows() == Matrix_right.num_rows())  &&  (this->num_columns() == Matrix_right.num_columns()) );
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = 1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasSaxpy(size_matrix,alpha, Matrix_right.p(),incx,this->p(),incy );
	
}


template<>
void matrix_gpu<double>::plus_assign(matrix_gpu<double> const& Matrix_right)
{
	assert( (this->num_rows() == Matrix_right.num_rows())  &&  (this->num_columns() == Matrix_right.num_columns()) );
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = 1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasDaxpy(size_matrix,alpha, Matrix_right.p(),incx,this->p(),incy );
	
}


template<>
const matrix_gpu<float> operator+ ( matrix_gpu<float>& Matrix_left, const matrix_gpu<float>& Matrix_right)
{
	assert( (Matrix_left.num_rows() == Matrix_right.num_rows())  &&  (Matrix_left.num_columns() == Matrix_right.num_columns()) );
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_columns = Matrix_left.num_columns();
		
	matrix_gpu<float> Result(num_rows,num_columns);
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = 1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasSaxpy(size_matrix,alpha, Matrix_right.p(),incx,Matrix_left.p(),incy );
	
	return Result;
};


template<>
const matrix_gpu<double> operator + ( matrix_gpu<double>& Matrix_left, const matrix_gpu<double>& Matrix_right)
{
	assert( (Matrix_left.num_rows() == Matrix_right.num_rows())  &&  (Matrix_left.num_columns() == Matrix_right.num_columns()) );
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_columns = Matrix_left.num_columns();
	
	matrix_gpu<double> Result(num_rows,num_columns);
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = 1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasDaxpy(size_matrix,alpha, Matrix_right.p(),incx,Matrix_left.p(),incy );
	
	return Result;
};



/*----------------------------- substraction -----------------------------------------*/


template<>
void matrix_gpu<float>::minus_assign(matrix_gpu<float> const& Matrix_right)
{
	size_type incx = 1;
	size_type incy = 1;
	float alpha = 1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasSaxpy(size_matrix,alpha, Matrix_right.p(),incx,this->p(),incy );
	
}


template<>
void matrix_gpu<double>::minus_assign(matrix_gpu<double> const& Matrix_right)
{
	size_type incx = 1;
	size_type incy = 1;
	float alpha = -1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasDaxpy(size_matrix,alpha, Matrix_right.p(),incx,this->p(),incy );
	
}

template<>
const matrix_gpu<float> operator- ( matrix_gpu<float>& Matrix_left, const matrix_gpu<float>& Matrix_right)
{
	assert( (Matrix_left.num_rows() == Matrix_right.num_rows())  &&  (Matrix_left.num_columns() == Matrix_right.num_columns()) );
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_columns = Matrix_left.num_columns();
	
	matrix_gpu<float> Result(num_rows,num_columns);
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = -1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasSaxpy(size_matrix,alpha, Matrix_right.p(),incx,Matrix_left.p(),incy );
	
	return Result;
};


template<>
const matrix_gpu<double> operator - ( matrix_gpu<double>& Matrix_left, const matrix_gpu<double>& Matrix_right)
{
	assert( (Matrix_left.num_rows() == Matrix_right.num_rows())  &&  (Matrix_left.num_columns() == Matrix_right.num_columns()) );
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_columns = Matrix_left.num_columns();
	
	matrix_gpu<double> Result(num_rows,num_columns);
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = -1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasDaxpy(size_matrix,alpha, Matrix_right.p(),incx,Matrix_left.p(),incy );
	
	return Result;
};

/*----------------------------- svd -----------------------------------------*/

template<>
void svd(matrix_gpu<float> & M, matrix_gpu<float> & U, matrix_gpu<float> & V, matrix_gpu<float> & S )
{

	culaStatus s;
	
	char jobu = 'A';
	char jobvt = 'A';
	
	size_type num_rows = M.num_rows();
	size_type num_columns = M.num_columns();
		
	s = culaDeviceSgesvd(jobu,jobvt,num_rows,num_columns,M.p(),M.ld(),S.p(),U.p(),U.ld(),V.p(),V.ld());
	if(s != culaNoError)
	{
		std::cout << culaGetErrorInfo() << std::endl;
		std::cout << culaGetStatusString(s) << std::endl;
	}
	
};

/* the double version is not inside the basic version 
template<>
void svd(matrix_gpu<double> & M, matrix_gpu<double> & U, matrix_gpu<double> & V)
{
	
	culaStatus s;
	
	
	char jobu = 'A';
	char jobvt = 'A';
	
	size_type num_rows = M.num_rows();
	size_type num_columns = M.num_columns();
	
	vector_gpu<double> S(std::min(num_rows,num_columns),0);
	
	s = culaDeviceDgesvd(jobu,jobvt,num_rows,num_columns,M.p(),M.ld(),S.p(),U.p(),U.ld(),V.p(),V.ld());
	if(s != culaNoError)
	{
		std::cout << culaGetErrorInfo() << std::endl;
		std::cout << culaGetStatusString(s) << std::endl;
	}
	
};
 */

/*----------------------------- qr -----------------------------------------*/

template<>
void qr(matrix_gpu<float> & M,  matrix_gpu<float> & Q , matrix_gpu<float> & R)
{
	
	size_type num_rows = M.num_rows();
	size_type num_columns = M.num_columns();
	vector_gpu<float> TAU(std::min(num_rows,num_columns),0);
	
	culaStatus s;	
	s = culaDeviceSgeqrf(num_rows, num_columns, M.p(), M.ld(), TAU.p() );
	if(s != culaNoError)
	{
		std::cout << culaGetErrorInfo() << std::endl;
		std::cout << culaGetStatusString(s) << std::endl;
	}

};

/*----------------------------- << -----------------------------------------*/	
	
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

