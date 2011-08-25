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

#define  MAXORDER vli::size_poly_vli::value
typedef   vli_gpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value>  Vli_gpu;
typedef   vli_cpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value>  Vli_cpu;

int main (int argc, char * const argv[]) 
{
    gpu::gpu_manager* GPU;
    GPU->instance();
    
    Vli_cpu acpu, ccpu;
    fill_random(acpu);
    Vli_cpu bcpu(acpu);

    Vli_gpu agpu(acpu);
    Vli_gpu bgpu(acpu);
    Vli_gpu cgpu;


    acpu.negate();
    agpu.negate();
/*
    bcpu.negate();
    bgpu.negate();

*/

    std::cout << agpu  << std::endl;
    std::cout << acpu  << std::endl;


    ccpu = bcpu*acpu;
    cgpu = bgpu*agpu;



    std::cout << cgpu  << std::endl;
    std::cout << ccpu  << std::endl;
/*

       
    ccpu = acpu * bcpu;   
    cgpu = agpu * bgpu;   
       
    std::cout << agpu << " " << bgpu << std::endl;
    std::cout << acpu << " " << bcpu << std::endl;
    
    std::cout << cgpu << std::endl;
    std::cout << ccpu << std::endl;

  */
    /*
    monomial<Vli> ma(a);
    polynomial_cpu<Vli, MAXORDER> pa;    
    fill_poly_random(pa);

    //Init GPU
    vli_gpu< Vli::value_type,Vli::size> a_gpu(a);
    monomial<vli_gpu< Vli::value_type,Vli::size> > magpu(a);
    polynomial_gpu<vli_gpu<Vli::value_type,Vli::size>, MAXORDER> pagpu(pa);

    pa = pa * ma; 
    pagpu = pagpu * magpu; 

    
    if(pagpu == pa){
       std::cout << "ok" << std::endl;
     }else{
        std::cout << "no ok" << std::endl;
     }
    */
    return 0;
}
