#include <iostream>
#include <cstdio>

#include <cassert>
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
typedef hp2c::monomial<large_int> monomial_type;
typedef hp2c::polynomial<large_int,11> polynomial_type;
typedef hp2c::polynomial<large_int,2*11> polynomial_typed;
typedef std::vector<polynomial_type> polynomial_vector_type;
#define SIZEV 16384
int main (int argc, char * const argv[]) 
{
    polynomial_vector_type v1gmp(SIZEV);
    polynomial_vector_type v2gmp(SIZEV);
    polynomial_type pgmp;
    polynomial_typed pgmpd;
    
#ifdef VLI_USE_GPU
    gpu::gpu_manager* gpu;
    gpu->instance();
#endif
    
    vector_type_cpu v1(SIZEV);
    vector_type_cpu v2(SIZEV);
    polynomial_result_type_cpu result_pure_cpu,result_mix_cpu_gpu,  result_cpu_gpu  ;
    
    fill_vector_random(v1,2);
    fill_vector_random(v2,2);

/*
     for (int i =0 ; i < SIZEV; ++i)
         for(int j = 0; j < 11; j++)
             for(int k = 0; k < 11; k++){
                 v1gmp[i](j,k) = v1[i](j,k).get_str();
                 v2gmp[i](j,k) = v2[i](j,k).get_str();
              }
*/
    TimerOMP t1("CPU openmp");
    t1.begin();
    result_pure_cpu = vli::detail::inner_product_openmp(v1,v2);
    t1.end();
/*
    TimerOMP t2("CPU gmp");
    t2.begin();
       pgmpd = inner_product(v1gmp,v2gmp);
    t2.end();
*/
#ifdef VLI_USE_GPU
    TimerOMP t3("MIX CPU/GPU openmp");
    t3.begin();    
    result_mix_cpu_gpu = vli::detail::inner_product_openmp_gpu(v1,v2);
    t3.end();
    
    if(result_mix_cpu_gpu ==result_pure_cpu ) {printf("OK \n"); } else{printf("NO OK \n"); }  
#endif
/*
    std::cout << result_pure_cpu << std::endl;
    std::cout << result_mix_cpu_gpu << std::endl;
*/
    return 0;
}


