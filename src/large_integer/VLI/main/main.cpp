#include <iostream>
#include <cstdio>
/*
#include "use_gmp_integers.hpp"
#include "minimal_polynomial.hpp"
*/
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


#define SIZE 21

using vli::vli_cpu;
using vli::max_int_value;
using vli::monomial;
using vli::polynomial_cpu;
using vli::vector_polynomial_cpu;

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::test::fill_vector_random;


typedef vli_cpu<long unsigned int, 3> vli_type_cpu;

typedef vli::monomial<vli_type_cpu> monomial_type_cpu;

typedef vli::polynomial_cpu< vli_type_cpu, SIZE > polynomial_type_cpu;
typedef vli::polynomial_cpu< vli_type_cpu, 2*SIZE > polynomial_result_type_cpu;

typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;
/*
typedef mpz_class large_int;
typedef hp2c::monomial<large_int> monomial_type;
typedef hp2c::polynomial<large_int> polynomial_type;
typedef std::vector<polynomial_type> polynomial_vector_type;
*/
int main (int argc, char * const argv[]) 
{
    
   
#ifdef VLI_USE_GPU
    gpu::gpu_manager* gpu;
    gpu->instance();
#endif
    /*
    vli_type_cpu a,b;
    vli::test::fill_random(a);
    vli::test::fill_random(b);
    
    b*=a;
    
    std::cout << b << std::endl;
    */
    vector_type_cpu v1(256);
    vector_type_cpu v2(256);
    polynomial_result_type_cpu result_pure_cpu,result_mix_cpu_gpu;
    
/*
    polynomial_vector_type v1gmp(1024);
    polynomial_vector_type v2gmp(1024);

    polynomial_result_type_cpu result;
    polynomial_result_type_cpu result_acc;
    polynomial_result_type_cpu result_plain;
    polynomial_result_type_cpu result_mix_cpu_gpu;
    polynomial_result_type_cpu result_pure_cpu; 
*/
    fill_vector_random(v1,1);
    fill_vector_random(v2,1);
/*
    // b as block
    polynomial_type_cpu p1,p2;
    polynomial_result_type_cpu r1,rb1;
    
    vli::test::fill_poly_random(p1);
    vli::test::fill_poly_random(p2);

    polynomial_type_cpu pb1(p1),pb2(p2);

    Timer A("classic");
    A.begin();
    poly_multiply(r1,p1,p2);
    A.end();

    Timer B("new");
    B.begin();
    poly_multiply_block(rb1,pb1,pb2);
    B.end();

    if(r1 == rb1) {std::cout << " ok " << std::endl;}
 */
 

    Timer t1("CPU openmp");
    t1.begin();
    result_pure_cpu = vli::detail::inner_product_openmp(v1,v2);
    t1.end();
    
        
   //std::cout << result_pure_cpu << std::endl;

     Timer t2("MIX CPU/GPU openmp");
    t2.begin();    
    result_mix_cpu_gpu = vli::detail::inner_product_openmp_gpu(v1,v2);
    t2.end();
    
 //  std::cout << result_mix_cpu_gpu << std::endl;
    if(result_mix_cpu_gpu ==result_pure_cpu ) {printf("OK \n"); } else{printf("NO OK \n"); }  
    
/*
    Timer t2("MIX CPU/GPU");
    t2.begin();
    result = inner_product(v1,v2);
    t2.end();

    std::cout << result_pure_cpu << std::endl;

    std::cout << result << std::endl;

    

    TimerOMP t3("oldsyte");
    t3.begin();
    result_acc=  vli::detail::inner_product_accp(v1,v2);
    t3.end();

    TimerOMP t4("plain");
    t4.begin();
    result_plain=  vli::detail::inner_product_plain(v1,v2);
    t4.end();

#ifdef VLI_USE_GPU
    TimerOMP t2("GPU");
    t2.begin();
    vli::detail::inner_product_gpu_booster<vector_type_cpu> gpu_boost(v1,v2,v1.size());
    result_pure_gpu = polynomial_type_cpu(gpu_boost);
    t2.end();
#endif //VLI_USE_GPU


    std::cout << result << std::endl << std::endl <<std::endl;
#ifdef VLI_USE_GPU
    std::cout << result_pure_gpu << std::endl << std::endl << std::endl;
#endif //VLI_USE_GPU
    std::cout << result_pure_cpu << std::endl << std::endl << std::endl;
   
    if(
       result == result_pure_cpu
 #ifdef VLI_USE_GPU
 && result == result_pure_gpu
 #endif //VLI_USE_GPU
      )
        printf( "OK \n"); 
    else{
        printf( "NOT OK \n"); 
    }
*/
    return 0;
}


