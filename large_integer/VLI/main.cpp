
//#include "vli_config.h"
#include <iostream>
#include <cstdio>

#define SIZE_BITS 256


#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "monome/monome.h"
#include "monome/polynome.h"
#include "monome/polynome_gpu.h"

#include "utils/timings.h"

typedef int TYPE; 
using vli::vli_cpu;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial;
using vli::polynomial_gpu;

int main (int argc, char * const argv[]) 
{
 	gpu::gpu_manager* GPU;
	GPU->instance();
    vli_cpu<int,8> a;
    a[0] = 252;
    a[1] = 245;
    a[2] = 44;
    a[3] = 97;
    a[4] = 106;
    a[5] = 219;
    a[6] = 198;
    a[7] = 0;
    
    vli_cpu<int,8> b = a+a+a;
    vli_cpu<int,8> c = a * vli_cpu<int,8>(3); 
	GPU->instance().destructor();
    return 0;
}

