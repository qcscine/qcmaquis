#include <iostream>
#include <cstdio>

#include "use_gmp_integers.hpp"
#include "minimal_polynomial.hpp"

#include <boost/lexical_cast.hpp>
#ifdef VLI_USE_GPU
#include "vli/utils/gpu_manager.h"
#include "vli/utils/gpu_manager.hpp"
#endif //VLI_USE_GPU
#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial_cpu.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/vli_traits.hpp"
#include "utils/timings.h"
#include "regression/vli_test.hpp"

#include "vli/detail/kernels_cpu_gpu.hpp"


#define SIZE 11
#define Size 3

using vli::vli_cpu;
using vli::max_int_value;
using vli::monomial;
using vli::polynomial_cpu;
using vli::vector_polynomial_cpu;

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::test::fill_vector_random;


typedef vli_cpu< unsigned long int, 3> vli_type_cpu;
typedef vli_cpu< unsigned long int, 6> vli_result_type_cpu;


typedef vli::monomial<vli_type_cpu> monomial_type_cpu;

typedef vli::polynomial_cpu< vli_type_cpu, SIZE > polynomial_type_cpu;
typedef vli::polynomial_cpu< vli_result_type_cpu, 2*SIZE > polynomial_result_type_cpu;

typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;

typedef mpz_class large_int;

int main (int argc, char * const argv[]) 
{
    /*
    polynomial_type_cpu p1,p2,pres;
    polynomial_result_type_cpu pres2;
    vli::test::fill_poly_random(p1);
    vli::test::fill_poly_random(p1);

    pres2 = p1 * p2;
    */
/*
    vli_type_cpu a,b;
    vli_result_type_cpu c;
    
    vli::test::fill_random(a,3);
    vli::test::fill_random(b,3);
    
    vli::detail::kernels_multiplication_classic<unsigned long int,Size>(&c[0],&a[0],&b[0]);

//    c = a * b;
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    mpz_class cgmp = agmp * bgmp;

    if(c.get_str() == cgmp.get_str()) {std::cout << " ok " << std::endl;}
    */
    
#ifdef VLI_USE_GPU
    gpu::gpu_manager* gpu;
    gpu->instance();
#endif
    
    vector_type_cpu v1(16384);
    vector_type_cpu v2(16384);
    polynomial_result_type_cpu result_pure_cpu,result_mix_cpu_gpu,  result_cpu_gpu  ;
    
    fill_vector_random(v1,1);
    fill_vector_random(v2,1);
    
    TimerOMP t1("CPU openmp");
    t1.begin();
    result_pure_cpu = vli::detail::inner_product_openmp(v1,v2);
    t1.end();

#ifdef VLI_USE_GPU
    TimerOMP t3("MIX CPU/GPU openmp");
    t3.begin();    
    result_mix_cpu_gpu = vli::detail::inner_product_openmp_gpu(v1,v2);
    t3.end();
    
    if(result_mix_cpu_gpu ==result_pure_cpu ) {printf("OK \n"); } else{printf("NO OK \n"); }  
#endif


    return 0;
}


