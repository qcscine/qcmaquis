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
#include "assert.h"
#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"

/* Notes */
/*To be change in futur = lda will be integrated in ETHZ matrix class */


namespace gpu
{

/**
* this class just initialized/shut down the cublas and cula library
*/	
class Simu
{
public:

/**
* the constructor starts cublas or cula, note : cula initialized also cublas ...
*/ 
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

/**
* The destructor shut down cublas or cula ...
*/ 	
	~Simu()
	{
#ifdef __CUBLAS__		
		cublasShutdown();		
#endif
		
#ifdef __CULA__		
		culaStatus s;
		culaShutdown();	
		if(s != culaNoError)
		{
			std::cout << " end : " << culaGetErrorInfo() << std::endl;
		}
#endif//
	}
};

/**
* GPU matrix class (pure GPU or mix mode CPU/GPU). It is based on the cublas library
*/
template<class T> 
class matrix_gpu
{	
public:
	
/**
* dummy constructor, no allocation, no initialization
*/
	matrix_gpu():size1_(0), size2_(0), ld_(0)
	{
	}
	
/**
* Constructor m(rows) X n(columns) matrix without given leading dimension, only use for the pure GPU blas operation.
* In this case, the leading dimension (ls) is equal to the number of row.
* We only allocates the memory on the GPU
*/
	matrix_gpu(size_type size1, size_type size2):size1_(size1), size2_(size2), ld_(size1) 	
	{
		cublasAlloc( size1*size2, sizeof(T), (void**)&p_ );
	}

	
/**
* Constructor m(rows) X n(columns) matrix with given leading dimension, only use for the mix CPU/GPU blas operation.
* In this case, a leading dimesion is given by the constructor argument.
* We only allocates the memory on the GPU
*/ 	
	matrix_gpu(size_type num_rows, size_type num_cols, size_type ld):size1_(num_rows), size2_(num_cols), ld_(ld) 	
	{
		cublasAlloc( num_rows*num_cols, sizeof(T), (void**)&p_ );
	}

/**
* Constructor m(rows) X n(columns) matrix without given leading dimension (ld = # of rows), 
* The memory is allocated on the GPU and initialized by a initial value.
*/ 	
	matrix_gpu(size_type size1, size_type size2, T value):size1_(size1), size2_(size2), ld_(size1) //To be change in futur	
	{
		size_type size_matrix = size1*size2;
		stat_ = cublasAlloc( size_matrix, sizeof(T), (void**)&p_ );		
		assert(true == CheckError(" cudaMalloc constructor matrix"));

		std::vector<T>  Array(size1*size2, value);		
		cublasSetMatrix(size1,size2,sizeof(T),&Array[0],size1,p(),size1);
		assert(true == CheckError(" cublasSetMatrix constructor matrix"));
	}
/**
* Copy constructor
*/
    matrix_gpu(matrix_gpu const& r)
        :size1_(r.size1_), size2_(r.size2_),ld_(r.ld_)
    {
        cublasAlloc(size1_*size2_,sizeof(T), (void**) &p_ );
        CheckError(" cudaMalloc copy constructor matrix");
		cudaMemcpy( p_, r.p_, size1_*size2_*sizeof(T) , cudaMemcpyDeviceToDevice);
    }
	
/**	
* Constructor m(rows) X n(columns) matrix using a blas::dense_matrix, m(# of rows), n(# of columns) and ld parameters. 
* These parameters can be different of the blas::dense_matrix, in case of divising the matrix B (from C = AXB)for
* the CPU/GPU mix mode
*/
	template<class MemoryBlock>
	matrix_gpu(blas::dense_matrix<T, MemoryBlock> const & Matrix_cpu, size_type nrows, size_type ncols, size_type ld)
	{
		size_type size_matrix = Matrix_cpu.num_columns()*Matrix_cpu.num_rows();
		stat_ = cublasAlloc( size_matrix, sizeof(T), (void**)&p_ );	
		assert(true == CheckError(" cudaMalloc constructor matrix"));
		cublasSetMatrix (nrows, ncols, sizeof(T),  &Matrix_cpu(0,0), Matrix_cpu.stride2(),p(), ld);	
	};
	
/**	
* Constructor m(rows) X n(columns) random matrix without given leading dimension, only for dev
*/
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
	
/**
* Copy constructor from CPU Matrix to GPU Matrix, use in pure GPU mode 
*/
	template<class MemoryBlock>
	matrix_gpu(blas::dense_matrix<T, MemoryBlock> const & Matrix_cpu):size1_(Matrix_cpu.num_rows()),size2_(Matrix_cpu.num_columns()),ld_(Matrix_cpu.stride2())
	{
		size_type size_matrix = size1_*size2_;
		stat_ = cublasAlloc( size_matrix, sizeof(T), (void**)&p_ );	
		assert(true == CheckError(" cudaMalloc constructor matrix"));
		
		cublasSetMatrix(size1_,size2_,sizeof(T), &Matrix_cpu(0,0),size1_,p_,ld_);
		assert(true == CheckError(" cublasSetMatrix constructor matrix"));
		
	};

	
/**
* Copy fonction from CPU Matrix to GPU Matrix, be carefull the memory on the GPU must be allocated before !
* If the GPU matrix exists, it will be overload.
*/	
	template <typename MemoryBlock>
	void copy_matrix_from_cpu(blas::dense_matrix<T,MemoryBlock> const& m_cpu)
	{
		assert( size1_ == num_rows(m_cpu) );
		assert( size2_ == num_columns(m_cpu) );
		cublasSetMatrix(size1_,size2_,sizeof(T),p(),ld_,&m_cpu(0,0),m_cpu.stride2());
	}
	
/**
* Copy fonction from GPU Matrix to CPU Matrix, be carefull the memory on the CPU must be allocated before !
* If the CPU matrix exists, it will be overload.
*/		
	template <typename MemoryBlock>
	void copy_matrix_to_cpu(blas::dense_matrix<T,MemoryBlock>& m_cpu) const
	{
		assert( size1_ == blas::num_rows(m_cpu));
		assert( size2_ == blas::num_columns(m_cpu) );
		cublasGetMatrix(size1_,size2_,sizeof(T),p(),ld_,&m_cpu(0,0),m_cpu.stride2());			
	}
	
/**
* overloading of the operator =, to modelize the CPU_Matrix = GPU_Matrix
*/
	template <typename MemoryBlock>
	operator blas::dense_matrix<T,MemoryBlock>()
	{
		blas::dense_matrix<T,MemoryBlock> m(size1_,size2_);
		copy_matrix_to_cpu(m);
		return m;
	}

/**
* destructor, we deallocate the GPU Memory 
*/
	~matrix_gpu()
	{
		cublasFree(p_);
	}

/**
* Overload the operator = it copies from GPU to GPU, be cautious of the size of your card
*/	
	matrix_gpu& operator = (matrix_gpu m)
	{
        swap(m);
		return *this;
	}

/**
* Swaps all content of two matrices
*/ 
    void swap(matrix_gpu& m)
    {
        std::swap(size1_,m.size1_);
        std::swap(size2_,m.size2_);
        std::swap(ld_,m.ld_);
        std::swap(p_,m.p_);
    }
/**
* Swaps all content of two matrices
*/ 
    friend void swap(matrix_gpu& m1, matrix_gpu& m2)
    {
        m1.swap(m2);
    }

/**
* Initialize a Identity Matrix, the GPU must be allocated before
*/
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

/**
* return a const pointer of the data inside the GPU memory, keep in mind, a reading/writing direct access is impossible  
*/
	inline const T* p () const
	{
		return p_;
	}

/**
* return a pointer of the data inside the GPU memory, keep in mind, a reading/writing direct access is impossible  
*/	
	inline  T* p () 
	{
		return p_;
	}

/**
* return a reference on the number of rows
*/
	size_type& size1() 
	{
		return size1_;
	}

/**
* return a reference on the number of columns
*/	
	size_type& size2() 
	{
		return size2_;
	}

/**
* return the number of rows, boost notation
*/	
	inline const size_type size1() const
	{
		return size1_;
	}

/**
* return the number of columns, boost notation
*/		
	inline const size_type size2() const
	{
		return size2_;
	}
	
/**
* return the number of rows, classical notation
*/		
	inline const size_type num_rows() const
	{
		return size1_;
	}

/**
* return the number of columns, classical notation
*/			
	inline const size_type num_columns() const
	{
		return size2_;
	}

/**
* return the leading dimension
*/			
	inline const size_type ld() const
	{
		return ld_;
	}

/**
* return a reference on the leading dimension
*/		
	size_type& ld() 
	{
		return ld_;
	}

/**
* overload the += operator
*/	
	matrix_gpu<T>& operator += (matrix_gpu<T> const& Matrix_right) 
	{
		plus_assign(*this,Matrix_right);
		return *this;
	}
	
/**
* overload the -= operator
*/	
	matrix_gpu<T>& operator -= (matrix_gpu<T> const& Matrix_right) 
	{
		minus_assign(*this,Matrix_right);
		return *this;
	}

/**
* overload the *= operator
*/	
	matrix_gpu<T>& operator *= (matrix_gpu<T> const& Matrix_right)
	{
		multiplies_assign(*this, Matrix_right);
		return *this;
	}

/** 
* Check the error of Cublas library, only works if #define NDEBUG is not defined in the main.
* It is build with the assert function
*/ 
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
/**
* Used inside the *= overload
*/
	void multiplies_assign (matrix_gpu<T> const& Matrix_right);

/**
* Used inside the += overload
*/	
	void plus_assign(matrix_gpu<T> const& Matrix_right);

/**
* Used inside the -= overload
*/
	void minus_assign(matrix_gpu<T> const& Matrix_right);
	
private:

/**
* pointer on the data inside the GPU, keep in mind a direct access to the datas
* is impossible 
*/
	T* p_;

/**	
* number of raws,
*/
	size_type size1_;

/**	
* number of columns
*/	
	size_type size2_;

/**	
* leading dimension
*/	
	size_type ld_;

/**
* cublas status
*/ 
	cublasStatus stat_;
};
	
template <class T>
void plus_assign(matrix_gpu<T>& Matrix_this, matrix_gpu<T> const& Matrix_right)
{
	Matrix_this.plus_assign(Matrix_right);
}

template <class T>
void minus_assign(matrix_gpu<T>& Matrix_this, matrix_gpu<T> const& Matrix_right)
{
	Matrix_this.minus_assign(Matrix_right);
}

template <class T>
void multiplies_assign(matrix_gpu<T>& Matrix_this, matrix_gpu<T> const& Matrix_right)
{
	Matrix_this.multiplies_assign(Matrix_right);
}
	
template <class T>
void multiplies_assign( matrix_gpu<T> &  Matrix_gpu, T const& t);
	
template <class T>
matrix_gpu<T> operator * ( matrix_gpu<T> const &  Matrix_left,  matrix_gpu<T>const & Matrix_right);	

// We keep exact notation of Andreas
template <class T>
void gemm( matrix_gpu<T>const & A, matrix_gpu<T> const & B, matrix_gpu<T>& C);	

/**
* matrix multiplication C = AxB, this function is used inside the hook function only for the pure GPU mode
*/
template <class T>
matrix_gpu<T> matrix_matrix_multiply( matrix_gpu<T>const & Matrix_left,  matrix_gpu<T>const & Matrix_right);
	
/**
* matrix multiplication C = AxB, this function is used inside the hook function only for the pure CPU/GPU mode
*/	
template <class T, class MemoryBlock>
void matrix_matrix_multiply(blas::dense_matrix<T,MemoryBlock>const & Matrix_left,  blas::dense_matrix<T,MemoryBlock>const & Matrix_right, blas::dense_matrix<T,MemoryBlock>& Matrix_result);
	
/**
* WARNING the Matrix_right will be overwrite for M_left + M_right
*/	
template<class T>
const matrix_gpu<T> operator + ( matrix_gpu<T>& Matrix_left, const matrix_gpu<T>& Matrix_right);

/**
* WARNING the Matrix_right will be overwrite for M_left - M_right  
*/	
template<class T>
const matrix_gpu<T> operator - ( matrix_gpu<T>& Matrix_left, const matrix_gpu<T>& Matrix_right);

template<class T>
void svd(matrix_gpu<T> & M, matrix_gpu<T> & U, matrix_gpu<T> & V, vector_gpu<T> & S);

template<class T>
void qr(matrix_gpu<T> & M,  matrix_gpu<T> & Q , matrix_gpu<T> & R);

template <class T>
std::ostream& operator<< (std::ostream& os, const  matrix_gpu<T> & Matrix_gpu);

	

}
#endif
