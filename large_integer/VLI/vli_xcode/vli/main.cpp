//
//  main.cpp
//  vli
//
//  Created by Tim Ewart on 15/06/11.
//  Copyright 2011 IBM. All rights reserved.
//

#include <iostream>
#define SIZE_BITS 256

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "engine/engine.h" 

genericbank* CPUBANK; 

int main (int argc, const char * argv[])
{
    CPUBANK = new cpubank();
    
    vli::vli_cpu<int> A;
    vli::vli_cpu<int> B;
    vli::vli_cpu<int> C;
    {
        vli::vli_cpu<int> D;
    }
    vli::vli_cpu<int> D;

    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

