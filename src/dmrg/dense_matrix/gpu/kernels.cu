
/**
Do not forget NVCC will compile this file, it dislikes STL and boost ....
Just very basic C++ or C (template ok). 
*/

/**
	The link between the DMRG code (CPU) and the DMRG code (GPU) is done inside the file  list_functions.h
	where the definition of function/wrapper are presented
*/
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
transpose kernel, to be tune
*/

/*
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
*/

template <typename T>
__global__ void transposeNaive(T odata, T idata, int num_rows, int num_columns, int ld)
{
   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   int index_in, index_out;
   
   __shared__ float block[16][16];
   
   
   	int qx=0;
	int rx=0;	
	qx = num_columns/blockDim.x;
	rx = num_columns - qx*blockDim.x ;
	
	int qy=0;
	int ry=0;	
	qy = num_rows/blockDim.y;
	ry = num_rows - qy*blockDim.y ;
	
   
   /*
	
   */
//	if (x < num_columns && y < num_rows) 	
//	{
	/*
		index_out = x + num_columns*y;
		index_in = y + num_rows*x;
	*/	
		index_in = x + num_columns*y;


		idata[y] = 999999 ;
		
		
//	}			

}



/**
swap row kernel
*/
template <typename T>
__global__ void SwapRows(T data, int num_rows, int num_columns, int ld, int i1, int i2)
{

	/**
	We optimize the swap with the local memory of the gpu card -> shared
	*/
   __shared__ float rows16[16];
   
   int niteration = 0;
   int nlimit = 0;
	/**
	The size of block must be smaller or equal to 16
	*/
   int dimensionofblocks = blockDim.x;
   
   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   
   /**
   local index
   */
   int index_i1, index_i2;

	/**
	Division Euclidian, need to adjust the transpose and avoid the overflow
	*/
   	int q=0;
	int r=0;	
	q = num_columns/blockDim.x;
	r = num_columns - q*blockDim.x ;

   if ((y == i1 || y == i2) && x < dimensionofblocks ) 	
   {
   
		nlimit = q*dimensionofblocks;
		for(niteration ; niteration < nlimit ; niteration += dimensionofblocks)
		{
			index_i2 = i2 + num_rows*(x+niteration);
			index_i1 = i1 + num_rows*(x+niteration); 

			rows16[x]   =  data[index_i1];
			data[index_i1] = data[index_i2];
			data[index_i2] = rows16[x];
		}
			
		/**
			second part of the swap done on the quotient r
			original test (r!=0 && x < r) to (x < r) (it should be enough)
		*/
		if(x < r)
		{
			index_i2 =  i2 + num_rows*(niteration+x);
			index_i1 =  i1 + num_rows*(niteration+x);

			rows16[x]      = data[index_i1];
			data[index_i1] = data[index_i2];
			data[index_i2] = rows16[x];
		}
   }	
}



/**
swap column kernel
*/
template <typename T>
__global__ void SwapColumns(T data, int num_rows, int num_columns, int ld, int i1, int i2)
{
	/**
	We optimize the swap with the local memory of the gpu card -> shared
	*/
   __shared__ float columns16[16];
   
   int niteration = 0;
   int nlimit = 0;
	/**
	The size of block must be smaller or equal to 16
	*/
   int dimensionofblocks = blockDim.y;
   
   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   
   /**
   local index
   */
   int index_i1, index_i2;

	/**
	Division Euclidian, need to adjust the swap and avoid the overflow
	*/
   	int q=0;
	int r=0;	
	q = num_rows/blockDim.y;
	r = num_rows - q*blockDim.y ;

   if ((x == i1 || x == i2) && y < dimensionofblocks ) 	
   {
		/**
			num_rows = q * dimensionofblocks + r
			first part of the swap done on the quotient q
			to do or not to do : unroll ?
		*/
		nlimit =q*dimensionofblocks;
		
		for(niteration ; niteration < nlimit ; niteration += dimensionofblocks)
		{
			index_i2 = i2*num_columns + y + niteration;
			index_i1 = i1*num_columns + y + niteration;

			columns16[y]   =  data[index_i1];
			data[index_i1] = data[index_i2];
			data[index_i2] = columns16[y];
		}
		
		/**
			second part of the swap done on the quotient r
			original test (r!=0 && y < r) to (y < r) (it should be enough)
		*/
		if(y < r)
		{
			index_i2 = i2*num_columns + y + niteration;
			index_i1 = i1*num_columns + y + niteration;

			columns16[y]   = data[index_i1];
			data[index_i1] = data[index_i2];
			data[index_i2] = columns16[y];
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

#ifdef NDEBUG
	printf("%i %i \n",dimgrid.x, dimgrid.y);
	printf("%i %i \n",dimblock.x, dimblock.y);
	printf("%i %i \n",dimthread.x, dimthread.y);
#endif
	SwapColumns <<< dimgrid, dimblock >>>(A_GPU_p, int(num_rows), int(num_columns), int(ld), int(i1), int(i2));
}



//void swap_columns(size_type j1, size_type j2);

