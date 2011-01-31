
/**
Do not forget NVCC will compile this file, it dislikes STL and boost ....
Just very basic C++ or C (template ok). 
*/

/**
	The link between the DMRG code (CPU) and the DMRG code (GPU) is done inside the file  list_functions.h
	where the definition of function/wrapper are presented
*/
#include "list_functions.h"
#include <iostream>

/**
the size of each block is 16X16
*/
#define NUM 16;

/**
	Do not forget, change to the desired type 
*/
typedef float TYPE;


/**
transpose kernel, to be tune
index_in and index_out is a mess => due to the row and columns order convention
*/
template <typename T>
__global__ void transposeNaive(T odata, T idata, int num_rows, int num_columns, int ld)
{
   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   int index_in, index_out;
   
	if (x < num_columns && y < num_rows) 	
	{
		index_out = x + num_columns*y;
		index_in = y + num_rows*x;
		odata[index_out] = idata[index_in] ;
	}			
}

/**
swap row kernel
*/
template <typename T>
__global__ void SwapNaiveRows(T data, int num_rows, int num_columns, int ld, int i1, int i2)
{

   __shared__ float rows[16];
   
   int numberofblocks = gridDim.x;

	/**
	The size of block must be smaller or equal to 16
	*/
   int dimensionofblocks = blockDim.x;
   
   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   int index_i1, index_i2;
   

   if ((x == i1 || x == i2) && y < dimensionofblocks) 	
   {
		for(int i=0;i<numberofblocks;i++)
		{
				
		index_i1 = i1 + num_rows*(y+i*dimensionofblocks);  // (y +i*dimensionofblocks) ;
		index_i2 = i2 + num_rows*(y+i*dimensionofblocks);
		
		rows[y]        =  data[index_i1];
		data[index_i1] = data[index_i2];
		data[index_i2] =  rows[y];

		}				
   }	
}



/**
swap column kernel
*/
template <typename T>
__global__ void SwapNaiveColumns(T data, int num_rows, int num_columns, int ld, int i1, int i2)
{

   __shared__ float rows[16];
   
   int numberofblocks = gridDim.y;

	/**
	The size of block must be smaller or equal to 16
	*/
   int dimensionofblocks = blockDim.y;
   
   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   int index_i1, index_i2;
   

   if ((x == i1 || x == i2) && y < dimensionofblocks ) 	
   {
		for(int i=0;i<numberofblocks;i++)
		{
			index_i2 = i2*num_columns + y + i*(dimensionofblocks);
			index_i1 = i1*num_columns + y + i*(dimensionofblocks);
			if(index_i1 < num_rows)
			{
				rows[y]        =  data[index_i1];
				data[index_i1] = data[index_i2];
				data[index_i2] = rows[y];
			}
		}
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
	transposeNaive <<< dimgrid, dimblock >>>(B_GPU_p, A_GPU_p, int(num_rows), int(num_columns), int(ld));

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
	SwapNaiveRows <<< dimgrid, dimblock >>>(A_GPU_p, int(num_rows), int(num_columns), int(ld), int(i1), int(i2));
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

#ifdef NDEBUG
	printf("%i %i \n",dimgrid.x, dimgrid.y);
	printf("%i %i \n",dimblock.x, dimblock.y);
	printf("%i %i \n",dimthread.x, dimthread.y);
#endif
	SwapNaiveColumns <<< dimgrid, dimblock >>>(A_GPU_p, int(num_rows), int(num_columns), int(ld), int(i1), int(i2));
}



//void swap_columns(size_type j1, size_type j2);

