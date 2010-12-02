/*
 *  matrix_gpu.h
 *  
 *
 *  Created by Tim Ewart and Alex Kosemkov on 26.11.10.
 *  Copyright 2010 University of Geneva. All rights reserved.
 *
 */


#ifndef __MATRIX_GPU__
#define __MATRIX_GPU__

typedef std::size_t             size_type;

#include "cuda.h"
#include "cuda_runtime_api.h" // need for cudaMemcpy, I do not know why !!!!
#include "cula.h"
#include "cublas.h"
#include "vector_gpu.h"
#include "allocator.h"

/* Notes */
/*To be change in futur = lda will be integrated in ETHZ matrix class */


namespace gpu
{

class Simu
{
public:
	
	Simu()
	{
#ifdef __CUBLAS__	
		cuInit(0);
		cublasInit();
#endif
		
#ifdef __CULA__	
		culaStatus s;
		culaInitialize();
		if(s != culaNoError)
		{
			std::cout << culaGetErrorInfo() << std::endl;
		}
#endif		
		
		
	}
	
	~Simu()
	{
#ifdef __CUBLAS__		
		cublasShutdown();		
#endif
		
#ifdef __CULA__		
		culaShutdown();	
#endif
	}
};


template<class T> 
class matrix_gpu
{	
public:
	
	//be cautious !!!! only use this constructor 

	matrix_gpu():size1_(0), size2_(0), ld_(0)
	{
	};
	
	matrix_gpu(size_type size1, size_type size2):size1_(size1), size2_(size2), ld_(size1) //To be change in futur	
	{
		cublasAlloc( size1*size2, sizeof(T), (void**)&p_ );
	};

	matrix_gpu(size_type size1, size_type size2, T value):size1_(size1), size2_(size2), ld_(size1) //To be change in futur	
	{
		size_type size_matrix = size1*size2;
		stat_ = cublasAlloc( size_matrix, sizeof(T), (void**)&p_ );		
		assert(true == CheckError(" cudaMalloc constructor matrix"));

		std::vector<T>  Array(size1*size2, value);		
		cublasSetMatrix(size1,size2,sizeof(T),&Array[0],size1,p(),size1);
		assert(true == CheckError(" cublasSetMatrix constructor matrix"));
	};

	//Only for me to generate a random matrix
	matrix_gpu(size_type size1, size_type size2, T value , bool boolean):size1_(size1), size2_(size2), ld_(size1) //To be change in futur	
	{
		size_type size_matrix = size1*size2;
		stat_ = cublasAlloc( size_matrix, sizeof(T), (void**)&p_ );		
		assert(true == CheckError(" cudaMalloc constructor matrix"));
		
		srand(value);
		std::vector<T ,alignment_allocator<T> >  Array(size1*size2);		
		
		for(int i = 0 ; i < size1*size2 ; i++)
			Array[i] = static_cast<T>(rand());
		
		cublasSetMatrix(size1,size2,sizeof(T),&Array[0],size1,p(),size1);
		assert(true == CheckError(" cublasSetMatrix constructor matrix"));
	};
	
	~matrix_gpu()
	{
		cublasFree(p_);
	}

	// Manual destructor could be usefull, just in case !
	void destroy()
	{
		this->~matrix_gpu();
	}
	

	matrix_gpu& operator = (const matrix_gpu& Matrix)
	{
		assert( (size1_ == Matrix.num_rows()) && (size2_ == Matrix.num_columns())  );
		
		size_type size_matrix = size1()*size2();
		
		cudaMemcpy(this->p_, Matrix.p(), size_matrix*sizeof(T) , cudaMemcpyDeviceToDevice);

		return *this;
	}

	
	void Identity()
	{
		assert( size1() == size2() );

		size_type size_matrix = size1()*size2();
		std::vector<T ,alignment_allocator<T> >  Array(size_matrix,0);
		
		for (int i=0; i < size1(); i++) 
			Array[i*size1()+i] = 1.0;
	
		cublasSetMatrix(size1(),size2(),sizeof(T),&Array[0],size1(),p(),size1());
		assert(true == CheckError(" cublasSetMatrix Identity"));
	
	}
	
	template<class M>
	void copy_cputogpu(M& Matrix_cpu)
	{
//		be cautious .... 
//		assert((size1_ == Matrix.num_rows() ) && (size2_ == Matrix.num_columns()));
		
		size_type size1 = Matrix_cpu.num_rows();
		size_type size2 = Matrix_cpu.num_columns();
	
		size1_ = size1 ; 
		size2_ = size2 ;
		ld_   = size1 ; //To be change in the futur
		
		cublasAlloc( size1*size2, sizeof(T), (void**)&p_ );
		cublasSetMatrix(size1,size2,sizeof(T), &Matrix_cpu(0,0),size1,p_,ld_);
	}

	
	inline const T* p () const
	{
		return p_;
	}
	
	inline  T* p () 
	{
		return p_;
	}

	
	//STL compatibility 
	inline const size_type size1() const
	{
		return size1_;
	}
	
	inline const size_type size2() const
	{
		return size2_;
	}
	
	inline const size_type num_rows() const
	{
		return size1_;
	}
	
	inline const size_type num_columns() const
	{
		return size2_;
	}
	
	inline const size_type ld() const
	{
		return ld_;
	}
	
	
	matrix_gpu<T>& operator += (matrix_gpu<T> const& Matrix_right) 
	{
		plus_assign(*this,Matrix_right);
		return *this;
	}
	
	matrix_gpu<T>& operator -= (matrix_gpu<T> const& Matrix_right) 
	{
		minus_assign(*this,Matrix_right);
		return *this;
	}
	
	matrix_gpu<T>& operator *= (matrix_gpu<T> const& Matrix_right)
	{
		multiplies_assign(*this, Matrix_right);
		return *this;
	}
	
	void plus_assign(matrix_gpu<T>& Matrix_this, matrix_gpu<T> const& Matrix_right)
	{
		Matrix_this.plus_assign(Matrix_right);
	}
	
	void minus_assign(matrix_gpu<T>& Matrix_this, matrix_gpu<T> const& Matrix_right)
	{
		Matrix_this.minus_assign(Matrix_right);
	}
	
	void multiplies_assign(matrix_gpu<T>& Matrix_this, matrix_gpu<T> const& Matrix_right)
	{
		Matrix_this.multiplies_assign(Matrix_right);
	}

	void multiplies_assign (matrix_gpu<T> const& Matrix_right);
	void plus_assign(matrix_gpu<T> const& Matrix_right);
	void minus_assign(matrix_gpu<T> const& Matrix_right);

	bool CheckError(std::string error)
	{
		switch (stat_) 
		{
			case CUBLAS_STATUS_NOT_INITIALIZED:
				std::cout << "CUBLAS_STATUS_NOT_INITIALIZED" + error << std::endl;
				return false;
				break;
			
			case CUBLAS_STATUS_MAPPING_ERROR:
				std::cout << "CUBLAS_STATUS_MAPPING_ERROR" + error << std::endl;
				return false;
				break;
	
			case CUBLAS_STATUS_INVALID_VALUE:
				std::cout << "CUBLAS_STATUS_INVALID_VALUE" + error << std::endl;
				return false;
				break;	
				
			default:
				std::cout << "CUBLAS_STATUS_SUCCESS" + error << std::endl;
				return true;
				break;
		}

	}
	
private:
	T* p_;
	//size1 and size2 to respect general_matrix implementation
	size_type size1_;
	size_type size2_;
	size_type ld_;
	cublasStatus stat_;
};

template <class T>
const matrix_gpu<T> operator * ( matrix_gpu<T>& Matrix_left, const matrix_gpu<T>& Matrix_right);		

/*  WARNING the Matrix_right will be overwrite for M_left + M_right or M_left - M_right  */
template<class T>
const matrix_gpu<T> operator + ( matrix_gpu<T>& Matrix_left, const matrix_gpu<T>& Matrix_right);

template<class T>
const matrix_gpu<T> operator - ( matrix_gpu<T>& Matrix_left, const matrix_gpu<T>& Matrix_right);
/* end WARNING */

template<class T>
void svd(matrix_gpu<T> & M, matrix_gpu<T> & U, matrix_gpu<T> & V, vector_gpu<T> & S);

template<class T>
void qr(matrix_gpu<T> & M,  matrix_gpu<T> & Q , matrix_gpu<T> & R);

template <class T>
std::ostream& operator<< (std::ostream& os, const  matrix_gpu<T> & Matrix_gpu);

	

}
#endif
