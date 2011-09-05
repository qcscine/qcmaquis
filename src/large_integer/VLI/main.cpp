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

typedef vli::monomial<vli_type_cpu> monomial_type_cpu;
typedef vli::monomial<vli_type_gpu> monomial_type_gpu;

typedef vli::polynomial_cpu< vli_type_cpu, vli::detail::size_poly_vli::value > polynomial_type_cpu;
typedef vli::polynomial_gpu< vli_type_gpu, vli::detail::size_poly_vli::value > polynomial_type_gpu;

typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;
typedef vli::vector_polynomial_gpu<polynomial_type_gpu> vector_type_gpu;


int main (int argc, char * const argv[]) 
{
    
    gpu::gpu_manager* GPU;
    GPU->instance();
       
    vector_type_cpu VaCPU( vli::detail::size_vector_vli::value); 
    
    vector_type_cpu result( vli::detail::size_vector_vli::value); 
    vector_type_gpu resultgpu( vli::detail::size_vector_vli::value); 

    fill_vector_random(VaCPU);
    vector_type_gpu VaGPU( VaCPU); 

    polynomial_type_cpu  A;
    fill_poly_random(A);

    polynomial_type_gpu  Ag(A);
  

    vli_type_cpu a;
    fill_random(a); 
    vli_type_gpu b(a);

    monomial_type_cpu mcpu(a,2,2);
    monomial_type_gpu mgpu(b,2,2);
    
    VaCPU[2] = VaCPU[3];//*monomial_type_cpu(1,0); 
    VaGPU[2] = VaGPU[3]; //*monomial_type_gpu(1,0); // VaGPU[2] = VaGPU[3] call 2 times the [] proxy operator
    VaCPU[1] = A;
    VaGPU[1] = Ag;
//    Ag = VaGPU[3];
    if(VaGPU == VaCPU ) {
       std::cout << "ok" << std::endl;
    }else{
       std::cout << "no ok" << std::endl;     
// yy      std::cout << Ag << std::endl;
    }

    

        
    //std::cout << VaCPU << std::endl;
//    std::cout << VaGPU << std::endl;

   
       
       
       
       
    return 0;
}


