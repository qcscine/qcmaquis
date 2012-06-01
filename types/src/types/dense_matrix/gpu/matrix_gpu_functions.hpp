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
//#include "cula.h"

#include "gpu/matrix_gpu.h"
#include "gpu/vector_gpu.h"

#include <utils/timings.h>

/*
My GT 330 does not support double so I develop, and debug with float.
run on CSCS with double, moreover we must respect the f77 philosophy.
*/


/*----------------------------- multiply -----------------------------------------*/
namespace gpu 
{
	
template <>
void multiplies_assign( matrix_gpu<float> &  Matrix_gpu, float const& t)
{
	size_type num = Matrix_gpu.num_rows()* Matrix_gpu.num_cols();
	cublasSscal(num, t , Matrix_gpu.p(), 1);
	cublasGetError();
}
	
template <>
void multiplies_assign( matrix_gpu<double> &  Matrix_gpu, double const& t)
{
	size_type num = Matrix_gpu.num_rows()* Matrix_gpu.num_cols();
	cublasDscal(num, t , Matrix_gpu.p(), 1);		
	cublasGetError();
}
	
	
template<>
void matrix_gpu<float>::multiplies_assign(matrix_gpu<float> const& Matrix_right)
{
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = this->num_rows();
	size_type num_cols = Matrix_right.num_cols();
	size_type K =  this->num_cols();
	
	cublasSgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_cols, K, 1, this->p(), this->ld(), Matrix_right.p(), Matrix_right.ld(), 0, this->p_, this->ld());	
	
}
	
template<>
void gemm(matrix_gpu<float> const& A, matrix_gpu<float> const& B, matrix_gpu<float> & C)
{
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = A.num_rows();
	size_type num_cols = B.num_cols();
	size_type K =  A.num_cols();
	cublasStatus s; 
	cublasSgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_cols, K, 1, A.p(), A.ld(), B.p(), B.ld(), 0, C.p(), C.ld());	
	if (s != CUBLAS_STATUS_SUCCESS )
	{
		maquis::cout << " Sgemm " <<cublasGetError() << std::endl;
	}

}
	
template<>
void matrix_gpu<double>::multiplies_assign(matrix_gpu<double> const& Matrix_right)
{
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = this->num_rows();
	size_type num_cols = Matrix_right.num_cols();
	size_type K =  this->num_cols();
	
	cublasDgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_cols, K, 1, this->p(), this->ld(), Matrix_right.p(), Matrix_right.ld(), 0, this->p_, this->ld());	
	
}

template<>
matrix_gpu<float> operator* ( matrix_gpu<float> const & Matrix_left,  matrix_gpu<float>const & Matrix_right)
{
	assert(Matrix_left.num_rows() == Matrix_right.num_cols());
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_cols = Matrix_right.num_cols();
	size_type K =  Matrix_left.num_cols();
	
	matrix_gpu<float> Result(num_rows,num_cols);
	
	cublasSgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_cols, K, 1, Matrix_left.p(), Matrix_left.ld(), Matrix_right.p(), Matrix_right.ld(), 0, Result.p(), Result.ld());
	
	return Result;
};

/**
* the mix mode CPU/GPU for SGEMM rule to respect : m_gpu = n_gpu + n_cpu 
* n_gpu = 8.n_cpu seems a good beginning to large matrix (superior to 10 000 ) 
*/
template<>
void matrix_matrix_multiply(maquis::types::dense_matrix<float,std::vector<float, std::allocator<float> > > const & lhs,maquis::types::dense_matrix<float,std::vector<float, std::allocator<float> > > const & rhs,maquis::types::dense_matrix<float,std::vector<float, std::allocator<float> > >  & result_cpu)
{
	
	const char TRANS_LEFT  = 'N';
	const char TRANS_RIGHT = 'N';
	
//	const fortran_int_t m_gpu = lhs.num_rows() ; old
	const fortran_int_t m_gpu = rhs.num_cols() ; 
	
	const fortran_int_t n_cpu = static_cast<fortran_int_t> (0) ;
//	const fortran_int_t n_cpu = static_cast<fortran_int_t> (2061) ;
	/**
		note : n_cpu = m_gpu = 0 -> full GPU
			   n_cpu = m_gpu     -> full CPU
			   n_cpu = m_gpu/n   -> CPU/GPU
	*/
	
	//const fortran_int_t n_cpu = static_cast<fortran_int_t> (m_gpu) ;
	const fortran_int_t n_gpu = m_gpu - n_cpu ;
	
	const fortran_int_t m = lhs.num_rows(); 
	const fortran_int_t k = lhs.num_cols() ; 
	const fortran_int_t k_gpu = rhs.num_rows() ; 
	
	const float alpha = 1;
	const float beta = 0;	
	
	const fortran_int_t lhs_stride2 =   lhs.stride2();
	const fortran_int_t rhs_stride2 =   rhs.stride2();
	const fortran_int_t result_cpu_stride2 =   result_cpu.stride2();		
	
	gpu::matrix_gpu<float> lhs_gpu(lhs,m,k,m); 
	gpu::matrix_gpu<float> rhs_gpu(rhs,k,n_gpu,k_gpu); 
	 
	gpu::matrix_gpu<float> result_gpu(static_cast<size_t> (m),static_cast<size_t> (n_gpu),static_cast<size_t> (m_gpu)); 

	Timer tps("inside") ;
	
	tps.begin();	
	
	cublasSgemm('n', 'n', m, n_gpu, k, 1.0, lhs_gpu.p(), m,rhs_gpu.p(), k, 0.0, result_gpu.p(), m); 

	tps.end();
	
	sgemm_(&TRANS_LEFT,&TRANS_RIGHT,&m, &n_cpu, &k, &alpha, &lhs(0,0), &lhs_stride2 ,&rhs(0,0)+rhs.stride2()*n_gpu, &rhs_stride2, &beta,&result_cpu(0,0)+result_cpu.stride2()*n_gpu, &result_cpu_stride2); 
	
	cublasGetMatrix (m, n_gpu, sizeof(float), result_gpu.p(), m, &result_cpu(0,0), result_cpu.stride2()); 

}

	
template<>
void matrix_matrix_multiply(maquis::types::dense_matrix<double,std::vector<double, std::allocator<double> > > const & lhs,maquis::types::dense_matrix<double,std::vector<double, std::allocator<double> > > const & rhs,maquis::types::dense_matrix<double,std::vector<double, std::allocator<double> > >  & result_cpu)
{
	const char TRANS_LEFT  = 'N';
	const char TRANS_RIGHT = 'N';
	
//	const fortran_int_t m_gpu = lhs.num_rows() ; 
	const fortran_int_t m_gpu = rhs.num_cols() ; 
	const fortran_int_t n_cpu = static_cast<fortran_int_t> (m_gpu/7) ;
	const fortran_int_t n_gpu = m_gpu - n_cpu ;
	
	const fortran_int_t m = lhs.num_rows(); 
	const fortran_int_t k = lhs.num_cols() ; 
	const fortran_int_t k_gpu = rhs.num_rows() ; 
	
	const double alpha = 1;
	const double beta = 0;	
	
	const fortran_int_t lhs_stride2 =   lhs.stride2();
	const fortran_int_t rhs_stride2 =   rhs.stride2();
	const fortran_int_t result_cpu_stride2 =   result_cpu.stride2();		
	
	gpu::matrix_gpu<double> lhs_gpu(lhs,m,k,m); 
	gpu::matrix_gpu<double> rhs_gpu(rhs,k,n_gpu,k_gpu); 
	gpu::matrix_gpu<double> result_gpu(static_cast<size_t> (m),static_cast<size_t> (n_gpu),static_cast<size_t> (m_gpu)); 
	
	cublasDgemm('n', 'n', m, n_gpu, k, 1.0, lhs_gpu.p(), m,rhs_gpu.p(), k, 0.0, result_gpu.p(), m); 
	dgemm_(&TRANS_LEFT,&TRANS_RIGHT,&m, &n_cpu, &k, &alpha, &lhs(0,0), &lhs_stride2 ,&rhs(0,0)+rhs.stride2()*n_gpu, &rhs_stride2, &beta,&result_cpu(0,0)+result_cpu.stride2()*n_gpu, &result_cpu_stride2); 
	cublasGetMatrix (m, n_gpu, sizeof(double), result_gpu.p(), m, &result_cpu(0,0), result_cpu.stride2()); 
    maquis::cout << result_cpu << std::endl;
}
	

/* it does not work, I do not why
template<class MemoryBlock>
void matrix_matrix_multiply(maquis::types::dense_matrix<float, MemoryBlock> const & Matrix_left,maquis::types::dense_matrix<float, MemoryBlock> const & Matrix_right,maquis::types::dense_matrix<float, MemoryBlock>  & Matrix_result)
{
		
};

template<class MemoryBlock>
void matrix_matrix_multiply(maquis::types::dense_matrix<double, MemoryBlock> const & Matrix_left,maquis::types::dense_matrix<double, MemoryBlock> const & Matrix_right,maquis::types::dense_matrix<double, MemoryBlock>  & Matrix_result)
{
		
};
*/	
template<>
matrix_gpu<float> matrix_matrix_multiply(matrix_gpu<float> const & Matrix_left,matrix_gpu<float> const & Matrix_right)
{
//	assert(Matrix_left.num_rows() == Matrix_right.num_cols());
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_cols = Matrix_right.num_cols();
	size_type K =  Matrix_left.num_cols();
	
	matrix_gpu<float> Result(num_rows,num_cols);

	cublasSgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_cols, K, 1, Matrix_left.p(), Matrix_left.ld(), Matrix_right.p(), Matrix_right.ld(), 0, Result.p(), Result.ld());
	cublasGetError();

	
	return Result;
}

	template<>
	matrix_gpu<double> matrix_matrix_multiply(matrix_gpu<double> const & Matrix_left,matrix_gpu<double> const & Matrix_right)
	{
		//	assert(Matrix_left.num_rows() == Matrix_right.num_cols());
		
		char  TRANS_LEFT  = 'N';
		char  TRANS_RIGHT = 'N';
		
		size_type num_rows = Matrix_left.num_rows();
		size_type num_cols = Matrix_right.num_cols();
		size_type K =  Matrix_left.num_cols();
		
		matrix_gpu<double> Result(num_rows,num_cols);
		
		cublasDgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_cols, K, 1, Matrix_left.p(), Matrix_left.ld(), Matrix_right.p(), Matrix_right.ld(), 0, Result.p(), Result.ld());
		cublasGetError();
		
		
		return Result;
	}
	
template<>
matrix_gpu<double> operator * ( matrix_gpu<double> const & Matrix_left,  matrix_gpu<double> const & Matrix_right)
{
	assert(Matrix_left.num_rows() == Matrix_right.num_cols());
	
	char  TRANS_LEFT  = 'N';
	char  TRANS_RIGHT = 'N';
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_cols = Matrix_right.num_cols();
	size_type K =  Matrix_left.num_cols();
	
	matrix_gpu<double> Result(num_rows,num_cols);
	
	cublasDgemm( TRANS_LEFT, TRANS_RIGHT, num_rows, num_cols, K, 1, Matrix_left.p(), Matrix_left.ld(), Matrix_right.p(), Matrix_right.ld(), 0, Result.p(), Result.ld());
	
	return Result;
};

/*----------------------------- addition -----------------------------------------*/

template<>
void matrix_gpu<float>::plus_assign(matrix_gpu<float> const& Matrix_right)
{
	assert( (this->num_rows() == Matrix_right.num_rows())  &&  (this->num_cols() == Matrix_right.num_cols()) );
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = 1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasSaxpy(size_matrix,alpha, Matrix_right.p(),incx,this->p(),incy );
	
}


template<>
void matrix_gpu<double>::plus_assign(matrix_gpu<double> const& Matrix_right)
{
	assert( (this->num_rows() == Matrix_right.num_rows())  &&  (this->num_cols() == Matrix_right.num_cols()) );
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = 1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasDaxpy(size_matrix,alpha, Matrix_right.p(),incx,this->p(),incy );
	
}


template<>
const matrix_gpu<float> operator+ ( matrix_gpu<float>& Matrix_left, const matrix_gpu<float>& Matrix_right)
{
	assert( (Matrix_left.num_rows() == Matrix_right.num_rows())  &&  (Matrix_left.num_cols() == Matrix_right.num_cols()) );
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_cols = Matrix_left.num_cols();
		
	matrix_gpu<float> Result(num_rows,num_cols);
	
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
	assert( (Matrix_left.num_rows() == Matrix_right.num_rows())  &&  (Matrix_left.num_cols() == Matrix_right.num_cols()) );
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_cols = Matrix_left.num_cols();
	
	matrix_gpu<double> Result(num_rows,num_cols);
	
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
	assert( (Matrix_left.num_rows() == Matrix_right.num_rows())  &&  (Matrix_left.num_cols() == Matrix_right.num_cols()) );
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_cols = Matrix_left.num_cols();
	
	matrix_gpu<float> Result(num_rows,num_cols);
	
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
	assert( (Matrix_left.num_rows() == Matrix_right.num_rows())  &&  (Matrix_left.num_cols() == Matrix_right.num_cols()) );
	
	size_type num_rows = Matrix_left.num_rows();
	size_type num_cols = Matrix_left.num_cols();
	
	matrix_gpu<double> Result(num_rows,num_cols);
	
	size_type incx = 1;
	size_type incy = 1;
	float alpha = -1.0;
	size_type size_matrix = Matrix_right.size1()*Matrix_right.size2();
	
	cublasDaxpy(size_matrix,alpha, Matrix_right.p(),incx,Matrix_left.p(),incy );
	
	return Result;
};

/*----------------------------- svd -----------------------------------------*/
/*
template<>
void svd(matrix_gpu<float> & M, matrix_gpu<float> & U, matrix_gpu<float> & V, vector_gpu<float> & S )
{

	culaStatus s;
	
	char jobu = 'A';
	char jobvt = 'A';
	
	size_type num_rows = M.num_rows();
	size_type num_cols = M.num_cols();
		
	s = culaDeviceSgesvd(jobu,jobvt,num_rows,num_cols,M.p(),M.ld(),S.p(),U.p(),U.ld(),V.p(),V.ld());
	if(s != culaNoError)
	{
		maquis::cout << culaGetErrorInfo() << std::endl;
		maquis::cout << culaGetStatusString(s) << std::endl;
	}
	
};
*/	
	/*

// the double version is not inside the basic version 
template<>
void svd(matrix_gpu<double> & M, matrix_gpu<double> & U, matrix_gpu<double> & V)
{
	
	culaStatus s;
	
	
	char jobu = 'A';
	char jobvt = 'A';
	
	size_type num_rows = M.num_rows();
	size_type num_cols = M.num_cols();
	
	vector_gpu<double> S(std::min(num_rows,num_cols),0);
	
	s = culaDeviceDgesvd(jobu,jobvt,num_rows,num_cols,M.p(),M.ld(),S.p(),U.p(),U.ld(),V.p(),V.ld());
	if(s != culaNoError)
	{
		maquis::cout << culaGetErrorInfo() << std::endl;
		maquis::cout << culaGetStatusString(s) << std::endl;
	}
	
};
 */

/*----------------------------- qr -----------------------------------------*/
/*
template<>
void qr(matrix_gpu<float> & M,  matrix_gpu<float> & Q , matrix_gpu<float> & R)
{
	
	size_type num_rows = M.num_rows();
	size_type num_cols = M.num_cols();
	vector_gpu<float> TAU(std::min(num_rows,num_cols),0);
	
	culaStatus s;	
	s = culaDeviceSgeqrf(num_rows, num_cols, M.p(), M.ld(), TAU.p() );
	if(s != culaNoError)
	{
		maquis::cout << culaGetErrorInfo() << std::endl;
		maquis::cout << culaGetStatusString(s) << std::endl;
	}

};
*/
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

