
/**
Do not forget NVCC will compile this file, it dislikes STL and boost ....
Just very basic C++ or C (template ok). 
*/

/**
	The link between the DMRG code (CPU) and the DMRG code (GPU) is done inside the file  list_functions.h
	where the definition of function/wrapper are presented
*/
#include <iostream>
#include "list_functions.h"

/**
the size of each block is 16X16
*/
#define NUM 16;

/**
	Do not forget, change to the desired type inside the kernels because I use special memory
*/
typedef float TYPE;


/**
	this kernel comes from Nvidia nevertheless they consider the ROW order and we used the columns order
	therefore I have adapted it, it is a big mess between x,y and num_columns and num_rows.
*/
template <typename T>
__global__ void transpose(T odata, T idata, int num_rows, int num_columns, int ld)
{
		unsigned int BLOCK_DIM = 16;
       __shared__ float block[16][16+1];

        // read the matrix tile into shared memory
        unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
        unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
        if((yIndex < num_rows) && (xIndex < num_columns))
        {
                unsigned int index_in = xIndex * num_rows + yIndex;
                block[threadIdx.y][threadIdx.x] = idata[index_in];
        }

        __syncthreads();

        // write the transposed matrix tile to global memory
        xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
        yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
        if((yIndex < num_columns) && (xIndex < num_rows))
        {
                unsigned int index_out = xIndex * num_columns + yIndex;
                odata[index_out] = block[threadIdx.x][threadIdx.y];
        }		
}


/**
swap row kernel
*/
template <typename T>
__global__ void SwapRows(T data, int num_rows, int num_columns, int ld, int i1, int i2)
{
		__shared__ float rows_i1[16];
		__shared__ float rows_i2[16];

		int xIndex = blockIdx.x*blockDim.x + threadIdx.x;
		int yIndex = blockIdx.y*blockDim.y + threadIdx.y;

		if ((yIndex == i1) && (xIndex < num_columns))
		{
				unsigned int index_in = xIndex * num_rows + yIndex;
				rows_i1[threadIdx.x] = data[index_in];
		} 

		if ((yIndex == i2) && (xIndex < num_columns))
		{
				unsigned int index_in = xIndex * num_rows + yIndex;
				rows_i2[threadIdx.x] = data[index_in];
		} 
		
		__syncthreads();
		
		if ((yIndex == i1) && (xIndex < num_columns))
		{
				unsigned int index_in = xIndex * num_rows + yIndex;
				data[index_in] = rows_i2[threadIdx.x];
		} 
		
		
		if ((yIndex == i2) && (xIndex < num_columns))
		{
				unsigned int index_in = xIndex * num_rows + yIndex;
				data[index_in] = rows_i1[threadIdx.x];
		} 
}

/**
swap columns kernel
*/
template <typename T>
__global__ void SwapColumns(T data, int num_rows, int num_columns, int ld, int i1, int i2)
{
		__shared__ float rows_i1[16];
		__shared__ float rows_i2[16];

        int xIndex = blockIdx.x *blockDim.x + threadIdx.x;
        int yIndex = blockIdx.y *blockDim.y + threadIdx.y;

		if ((xIndex == i1) && (yIndex < num_rows))
		{
				unsigned int index_in = xIndex * num_rows + yIndex;
				rows_i1[threadIdx.y] = data[index_in];
		} 

		if ((xIndex == i2) && (yIndex < num_rows))
		{
				unsigned int index_in = xIndex * num_rows + yIndex;
				rows_i2[threadIdx.y] = data[index_in];
		} 
		
		__syncthreads();
		
		if ((xIndex == i1) && (yIndex < num_rows))
		{
				unsigned int index_in = xIndex * num_rows + yIndex;
				data[index_in] = rows_i2[threadIdx.y];
		} 
		
		
		if ((xIndex == i2) && (yIndex < num_rows))
		{
				unsigned int index_in = xIndex * num_rows + yIndex;
				data[index_in] = rows_i1[threadIdx.y];
		} 

}




/**
	Note that are the arguments float * A_GPU and float * B_GPU
	are the pointer _p of the class matrix_gpu
*/
template<>
void transpose(TYPE*  B_GPU_p, TYPE*  A_GPU_p, std::size_t num_rows, std::size_t num_columns, std::size_t ld)
{
	dim3 dimgrid;
	dim3 dimblock;
	dim3 dimthread;

	dimblock.x = NUM;
	dimblock.y = NUM;
	dimblock.z = 1;
		
	dimthread.x = NUM;
	dimthread.y = NUM;
	dimthread.z = 1;		
	
	dimgrid.x = (int(num_columns) + dimblock.x - 1)/ dimblock.x;
	dimgrid.y = (int(ld) + dimblock.y - 1)/ dimblock.y;
	dimgrid.z = 1;

#ifdef NDEBUG
	printf("%i %i \n",dimgrid.x, dimgrid.y);
	printf("%i %i \n",dimblock.x, dimblock.y);
	printf("%i %i \n",dimthread.x, dimthread.y);
#endif
	transpose <<< dimgrid, dimblock >>>(B_GPU_p, A_GPU_p, int(num_rows), int(num_columns), int(ld));

}

template<>
void swap_rows(TYPE* A_GPU_p, std::size_t num_rows, std::size_t num_columns, std::size_t ld , std::size_t i1, std::size_t i2)
{
	dim3 dimgrid;
	dim3 dimblock;
	dim3 dimthread;

	dimblock.x = NUM;
	dimblock.y = NUM;
	dimblock.z = 1;
		
	dimthread.x = NUM;
	dimthread.y = NUM;
	dimthread.z = 1;		
	
	dimgrid.x = (int(num_columns) + dimblock.x - 1)/ dimblock.x;
	dimgrid.y = (int(ld) + dimblock.y - 1)/ dimblock.y;
	dimgrid.z = 1;

#ifdef NDEBUG
	printf("%i %i \n",dimgrid.x, dimgrid.y);
	printf("%i %i \n",dimblock.x, dimblock.y);
	printf("%i %i \n",dimthread.x, dimthread.y);
#endif
	SwapRows <<< dimgrid, dimblock >>>(A_GPU_p, int(num_rows), int(num_columns), int(ld), int(i1), int(i2));
}


template<>
void swap_columns(TYPE* A_GPU_p, std::size_t num_rows, std::size_t num_columns, std::size_t ld , std::size_t i1, std::size_t i2)
{
	dim3 dimgrid;
	dim3 dimblock;
	dim3 dimthread;

	dimblock.x = NUM;
	dimblock.y = NUM;
	dimblock.z = 1;
		
	dimthread.x = NUM;
	dimthread.y = NUM;
	dimthread.z = 1;		
	
	dimgrid.x = (int(num_columns) + dimblock.x - 1)/ dimblock.x;
	dimgrid.y = (int(ld) + dimblock.y - 1)/ dimblock.y;	
	dimgrid.z = 1;	

//#ifdef NDEBUG
	printf("%i %i \n",dimgrid.x, dimgrid.y);
	printf("%i %i \n",dimblock.x, dimblock.y);
	printf("%i %i \n",dimthread.x, dimthread.y);
//#endif
	SwapColumns <<< dimgrid, dimblock >>>(A_GPU_p, int(num_rows), int(num_columns), int(ld), int(i1), int(i2));
}



//void swap_columns(size_type j1, size_type j2);

