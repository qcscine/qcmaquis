#include <iostream>
#include <cstdio>

#include <boost/lexical_cast.hpp>

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"

#include "polynomial/vector_polynomial_cpu.hpp"
#include "polynomial/vector_polynomial_gpu.hpp"
#include "polynomial/polynomial_gpu.hpp"
#include "polynomial/polynomial_cpu.hpp"
#include "polynomial/monomial.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "utils/timings.h"
#include "regression/common_test_functions.hpp"
#include "detail/vli_size_param.hpp"


using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial_cpu;
using vli::polynomial_gpu;
using vli::vector_polynomial_gpu;
using vli::vector_polynomial_cpu;

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::test::fill_vector_random;


typedef vli_cpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value> vli_type_cpu;
typedef vli_gpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value> vli_type_gpu;

typedef vli::monomial<vli_cpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value> > monomial_type_cpu;
typedef vli::monomial<vli_gpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value> > monomial_type_gpu;

typedef vli::polynomial_cpu<vli_cpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value>, vli::detail::size_poly_vli::value > polynomial_type_cpu;
typedef vli::polynomial_gpu<vli_gpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value>, vli::detail::size_poly_vli::value > polynomial_type_gpu;

typedef vli::vector_polynomial_gpu<polynomial_type_gpu> vector_type_gpu;
typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;


int main (int argc, char * const argv[]) 
{
    
    gpu::gpu_manager* GPU;
    GPU->instance();
    polynomial_type_cpu pcCPU;
    polynomial_type_gpu pcGPU;
    

    vector_type_cpu VaCPU(512),VbCPU(512); 
    
    fill_vector_random(VaCPU);
    fill_vector_random(VbCPU);

    vector_type_gpu VaGPU(VaCPU); 
    vector_type_gpu VbGPU(VbCPU); 


    Timer A("CPU");
    A.begin();            
    pcCPU =  inner_product(VaCPU,VbCPU);
    A.end();

    Timer B("GPU");
    B.begin();            
    pcGPU =  inner_product(VaGPU,VbGPU);
    B.end();
  
    if(pcGPU == pcCPU)
        std::cout << " OK " << std::endl;
    else{
        std::cout << pcGPU << std::endl;
        std::cout << pcCPU << std::endl;
    }
      
    return 0;
}


