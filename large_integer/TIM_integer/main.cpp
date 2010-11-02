#include <iostream>
#include <vector>
#include <emmintrin.h>
#include <OpenCL/OpenCL.h>

#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"

#include "CIntSIMD.h"
#include "CIntGPU.h"
#include "CintVL.h"



typedef  Clint<CContainerSIMD<std::size_t>, size_t> int128_SIMD;
typedef  Clint<CContainerGPU<std::size_t>, size_t> int128_GPU;

int main (int argc, char * const argv[])
{
	int128_SIMD A,B;
	int128_SIMD C(A);

	srand(3);
		
	A[0] = rand()%BASE;
	A[1] = rand()%BASE;

	C = A ;
	
	C = A + B;
	
    return 0;
}
