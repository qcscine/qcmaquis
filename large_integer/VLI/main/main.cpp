

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
#include "vli/polynomial/polynomial.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.h"
#include "vli/vli_traits.hpp"
#include "utils/timings.h"
#include "regression/vli_test.hpp"

#include "vli/detail/bit_masks.hpp"

#define Size1 3
#define Size2 6
#define Order 11

using vli::vli_cpu;
using vli::max_int_value;
using vli::monomial;
using vli::polynomial;
using vli::vector_polynomial;

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::test::fill_vector_random;
using vli::test::fill_vector_negate;


typedef vli_cpu< unsigned long int, Size1> vli_type_cpu;
typedef vli_cpu< unsigned long int, Size2> vli_result_type_cpu;


typedef vli::monomial<vli_type_cpu> monomial_type_cpu;

typedef vli::polynomial< vli_type_cpu, Order > polynomial_type_cpu;
typedef vli::polynomial< vli_result_type_cpu, 2*Order > polynomial_result_type_cpu;

typedef vli::vector_polynomial<polynomial_type_cpu> vector_type_cpu;

typedef mpz_class large_int;
typedef hp2c::monomial<large_int> monomial_type;
typedef hp2c::polynomial<large_int,Order> polynomial_type;
typedef hp2c::polynomial<large_int,2*Order> polynomial_typed;
typedef std::vector<polynomial_type> polynomial_vector_type;

template <typename VpolyVLI, typename VpolyGMP>
void InitPolyVLItoPolyGMP(VpolyVLI const& VVLI, VpolyGMP & VGMP)
{
    #pragma omp parallel for
    for (int i =0 ; i < (int)VVLI.size() ; ++i)
        for(int j = 0; j < Order; j++)
            for(int k = 0; k < Order; k++){
                VGMP[i](j,k) = VVLI[i](j,k).get_str();
            }
}

template <typename PolyVLI, typename PolyGMP>
bool ValidatePolyVLI_PolyGMP(PolyVLI const& PVLI, PolyGMP const& PGMP)
{
    bool b(true);
  #pragma omp parallel for
    for(std::size_t j = 0; j < PolyVLI::max_order; j++)
        for(std::size_t k = 0; k < PolyVLI::max_order; k++){
            if( PGMP(j,k).get_str() != PVLI(j,k).get_str()){
                 b = false;
            }
        }   
    return b;
}

int main (int argc, char * const argv[]) 
{
    
    int SizeVector =  atoi(argv[1]); 
     
 polynomial_vector_type v1gmp(SizeVector);
 polynomial_vector_type v2gmp(SizeVector);
 polynomial_type pgmp;
 polynomial_typed pgmpd;

#ifdef VLI_USE_GPU
 gpu::gpu_manager* gpu;
 gpu->instance();
 #endif

 vector_type_cpu v1(SizeVector);
 vector_type_cpu v2(SizeVector);

 polynomial_result_type_cpu result_pure_cpu  ;
 polynomial_result_type_cpu result_pure_cpu_omp  ;
 polynomial_result_type_cpu result_mix_cpu_gpu  ;
 polynomial_result_type_cpu result_cpu_gpu  ;
   
    fill_vector_random(v1,3);
    fill_vector_random(v2,2);

    fill_vector_negate(v1,3);
    fill_vector_negate(v2,2);

    InitPolyVLItoPolyGMP(v1,v1gmp);
    InitPolyVLItoPolyGMP(v2,v2gmp);

    Timer t3("CPU vli_omp");
    t3.begin();
      result_pure_cpu = vli::detail::inner_product_openmp(v1,v2);
    t3.end();

#ifdef VLI_USE_GPU
std::cout << " --------------------------- " << std::endl;
    TimerOMP t("GPU omp 1" );
    t.begin();   
      result_pure_cpu_omp = vli::detail::inner_product_gpu_omp(v1,v2);
    t.end();

std::cout << " --------------------------- " << std::endl;
    TimerOMP t1("GPU omp 2 ");
    t1.begin();   
      result_pure_cpu_omp = vli::detail::inner_product_gpu_omp(v1,v2);
    t1.end();
std::cout << " --------------------------- " << std::endl;
    TimerOMP t2("GPU omp 3 ");
    t2.begin();   
      result_pure_cpu_omp = vli::detail::inner_product_gpu_omp(v1,v2);
    t2.end();


    Timer t4("CPU gmp_omp");
    t4.begin();
       pgmpd = inner_product(v1gmp,v2gmp);
    t4.end();
#endif

/*
     std::cout << " GPU----------------------------------------------------- " << std::endl;
     std::cout << result_pure_cpu_omp << std::endl;
     std::cout << " CPU----------------------------------------------------- " << std::endl;
     std::cout << result_pure_cpu << std::endl;
*/
    if(result_pure_cpu_omp ==result_pure_cpu ) {printf("OK gpu\n"); } else{printf("NO OK gpu \n"); } 

    if(ValidatePolyVLI_PolyGMP(result_pure_cpu,pgmpd))
        std::cout << "validation GMP OK " << std::endl;

    return 0;
}



