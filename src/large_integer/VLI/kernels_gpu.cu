
/**
Do not forget NVCC will compile this file, it dislikes STL and boost ....
Just very basic C++ or C (template ok). 
*/

/**
	The link between the DMRG code (CPU) and the DMRG code (GPU) is done inside the file  list_functions.h
	where the definition of function/wrapper are presented
*/
#include <iostream>
#include "definition.h"

/**
the size of each block is 16X16
*/

/**
	Do not forget, change to the desired type inside the kernels because I use special memory
*/
typedef int TYPE;

/**
	C = A + B
	if this kernel we introduce new arithmetic
	A. Avizienis. Signed-digit number representations for fast parallel arithmetic. IRE Transactions on electronic computers, vol 10, pages 389-400, 1961
	No carry bit ^^'
	num_integer = # of vli
	ld = size of one vli
*/

extern const int NUM_SHARED = 128; // must be equal to the size of a block

template <typename T>
__global__ void addition_Avizienis_kernel_gpu(T x, T y , T z,int num_integer, int ld)
{

/**
	sharred array for the local calculation t must be larger (+1) due to the algo
*/

	__shared__ TYPE t[NUM_SHARED+1]; 
	__shared__ TYPE w[NUM_SHARED];
	__shared__ TYPE s[NUM_SHARED];

	int nsize = ld*num_integer; //remove if the size of the grid is exact 
	
/**
	just to remember how to calculate the index of the grid
	int xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
	int yIndex = blockIdx.y*blockDim.y + threadIdx.y; // all index on y
	int i = xIndex + yIndex* blockDim.x * gridDim.x ; // all index on xy 
*/
	int i = threadIdx.y + blockDim.y*threadIdx.x; // index on the block (0,0)
	int j = i; // the copy for the local memory


/**
	addition block by block to have a better access W/R of the shared memory
	We should have only one block on y; because y represents the size ot the vli
	
	k ( on x direction)
	--->
	______________________ ......
	|blockId|blockId|
	|(0,0)  |(1,0)	|
	|		|		|
	|		|		|
	|		|		|
	______________________ ......
	
	we make the addition on a block and we go to another block
	i is the index of the global memory
	j is the index of the shared memory
	
*/
	
/**
	inspire froan scan pattern design 
	http://developer.nvidia.com/object/cuda_training.html
	class : Stanford university CS193G Parallel Patterns I, slide 40
*/
	for(int k =0 ; k < gridDim.x ; k++) //loop on the block
	{
		if(i < nsize) // To avoid the corruption of the memory card (overflow), remove is the size of the grid is perfect
		{
			s[j] = 0;
			w[j] = 0;
			t[j] = 0;
		
			s[j] = x[i] + y[i];
		
			__syncthreads(); // wait on read

			// To do : optimize 
			if(s[j] > BASE_MINUS2)
				t[j+1] = 1;
			if(s[j] < MINUS_BASE_PLUS2)
				t[j+1] = -1;
			if(s[j]<=BASE_MINUS2 && s[j]>= MINUS_BASE_PLUS2)
				t[j+1] = 0;
			
			w[j] = s[j] - BASE*t[j+1];
	
			z[i] = w[j] + t[j];
		
			__syncthreads(); // wait on write
			
			i+=NUM_SHARED; // iterate to go to the next block
		}
	}
}

/**
	num_integer is the number of integer to construct the very large integer
	ld is the number of very large integer (size of vector)
*/

/**
classical addition, the operation " addition " is serial but we sum all number together 
*/
template <typename T>
__device__ void addition_kernel_gpu(T* x,   const T * y, T k)
{
	T carry_bit; 
	carry_bit = 0;
	*(x+k)    += *(y+k);
	carry_bit  = *(x+k) >> LOG_BASE;
	*(x+k)    %= BASE;
	*(x+k+1)  += carry_bit; //MAYBE PB TO CHECK
}

template <typename T>
__global__ void addition_classic_kernel_gpu(T* x, const T * y ,T num_integer, T ld)
{
	T xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
	T j = xIndex*ld; // index to be on the beginning of the vli (beginning of every columns)
	T k = 0;
	if(xIndex < num_integer) // the classical condition to avoid overflow
	{	
		for (size_type i = 0 ; i < ld; i++) 
		{
			k = i+j ;
			addition_kernel_gpu(x,y,k);
		}
	}

}

/**
classical multiplication, the operation " addition " is serial but we sum all number together 
*/

template <typename T>
__device__ void addition_kernel_gpu_noiter(T* x, const T* y )
{
	int carry_bit; 
	carry_bit = 0;
	*(x)    +=  *(y);
	carry_bit  = *(x) >> LOG_BASE;
	*(x)    %= BASE;
	*(x+1)  += carry_bit; //MAYBE PB TO CHECK
}

template <typename T>
__device__ void multiplication_kernel_up_gpu(const T*   x, const T*    y, T * r)	
{
	*r	    = ((*x & MASK_UP) >> LOG_BASE_HALF ) * (*y & MASK_DOWN);	
	*(r+1)	= ((*x & MASK_UP) >> LOG_BASE_HALF ) * ((*y & MASK_UP) >> LOG_BASE_HALF);
}		
	
template <typename T>
__device__ void multiplication_kernel_down_gpu(const T*  x,  const T*   y, T * r)	
{	
	*r     = (*x & MASK_DOWN) * (*y & MASK_DOWN);
	*(r+1) = (*x & MASK_DOWN) * ((*y & MASK_UP) >> LOG_BASE_HALF);
}

template <typename T>
__device__ void multiplication_kernel_base_reshaping_gpu(T*  a, const T*   b, T * r)	
{	
	int q1,q2;
	int r1,r2;
	q1 = q2 = r1 =r2 = 0;

	q1 = (*(a+1) + *b)/BASE_HALF;
	r1 = (*(a+1) + *b)%BASE_HALF;
	r1 = r1 * BASE_HALF;
	q2 = (r1 + *a)/BASE; 
	r2 = (r1 + *a)%BASE;
	*r  = r2;
	*(r+1) = q1 + q2 + *(b+1);
}

template <typename T>
__device__  void multiplication_block_gpu( const T*   x,   const  T*   y, T *r)	
{
	T a[2] = {0,0};
	T b[2] = {0,0};
	/**
	 Divide and conquer algo (see my notes - for the euclidian division tips)
		X <=> Xl Xr (half of the binary number)
	 x  Y <=> Yl Yr (half of the binary number)
	-------
	= 2^n XlYl + 2^(n/2) (XlYr + XrYl) + XrYr (multiplication_kernel_gpu_down and multiplication_kernel_gpu_up)
	------- 
	= (q1+q2 + Xl*Yl)*BASE + r2  (multiplication_kernel_base_reshaping)
	*/
	multiplication_kernel_down_gpu(x,y, a);
	multiplication_kernel_up_gpu(x,y, b);
	multiplication_kernel_base_reshaping_gpu(a,b, r);
	//r[0] = 9999;
}	

template <typename T>
__global__ void multiplication_classic_kernel_gpu(const T* x,  const T* y , T* z , int num_integer, int ld)
{
	T r[2] = {0,0};	//for local block calculation
	
	T xIndex = blockIdx.x*blockDim.x + threadIdx.x; // all index on x
	T i_ld = xIndex*ld; // index to be on the beginning of the vli (beginning of every columns)
	T m=0;
	
	if(xIndex < num_integer) // the classical condition to avoid overflow
	{
		//One of this two loops could be remove to do
		for (T i = 0 ; i < 1; i++) 
		{
			for(T j = 0 ; j < 1 ; j++)  			
			{	
				m = j + i;
				multiplication_block_gpu((x+i_ld+i), (y+i_ld+j), r);
				addition_kernel_gpu_noiter((z+i_ld+m),r);//,i_ld+m);	
				addition_kernel_gpu_noiter((z+i_ld+m+1),r+1);//,i_ld+m+1);					
				__syncthreads(); // Need, why (carry bit propagation) ?
			}
		}
	}
}


/*------------------------ MANAGEMENT ------------------------------------------------------------ */


void DeterminationGrid(dim3& dimgrid, dim3& dimblock, dim3& dimthread, int num_integer, int ld)
{
	dimblock.x = NUM;
	dimblock.y = ld;
	dimblock.z = 1;
		
	dimthread.x = NUM;
	dimthread.y = ld;
	dimthread.z = 1;		
	
	dimgrid.x = (int(num_integer) + dimblock.x - 1)/ dimblock.x;
	dimgrid.y = (int(ld) + dimblock.y - 1)/ dimblock.y;
	dimgrid.z = 1;
}

template<>
void addition_gpu(TYPE*  A,  const TYPE*  B, TYPE num_integer, TYPE ld)
{

	dim3 dimgrid;
	dim3 dimblock;
	dim3 dimthread;

	DeterminationGrid(dimgrid, dimblock, dimthread,num_integer,ld );
	
#ifdef NDEBUG
	printf("%i %i \n",dimgrid.x, dimgrid.y);
	printf("%i %i \n",dimblock.x, dimblock.y);
	printf("%i %i \n",dimthread.x, dimthread.y);
#endif

//	addition_Avizienis_kernel_gpu <<< dimgrid, dimblock >>>(A, B, C, num_integer, ld);
	addition_classic_kernel_gpu <<< dimgrid, dimblock >>>(A, B, num_integer, ld);
//	printf("%i", A);
}

template<>
void multiply_gpu(const TYPE*  A, const TYPE*  B, TYPE* C ,TYPE num_integer, TYPE ld)
{

	dim3 dimgrid;
	dim3 dimblock;
	dim3 dimthread;

	DeterminationGrid(dimgrid, dimblock, dimthread,num_integer,ld );
	
#ifdef NDEBUG
	printf("%i %i \n",dimgrid.x, dimgrid.y);
	printf("%i %i \n",dimblock.x, dimblock.y);
	printf("%i %i \n",dimthread.x, dimthread.y);
#endif

	multiplication_classic_kernel_gpu <<< dimgrid, dimblock >>>(A, B ,C, num_integer, ld);
}




