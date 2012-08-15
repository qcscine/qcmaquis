

#include <iostream>
#include <cstdio>
#include <cassert>


#include <boost/mpl/for_each.hpp>

#include "use_gmp_integers.hpp"
#include "minimal_polynomial.hpp"

#include <boost/lexical_cast.hpp>
#include "boost/tuple/tuple.hpp"

#ifdef VLI_USE_GPU
#include "vli/detail/gpu/inner_product_gpu_booster.hpp"
#include "vli/utils/gpu_manager.h"
#include "vli/utils/gpu_manager.hpp"
#endif //VLI_USE_GPU
#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial.hpp"
#include "vli/polynomial/variable.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/polynomial/polynomial_traits.hpp"
#include "vli/vli_cpu.h"
#include "vli/vli_traits.hpp"
#include "utils/timings.h"
//#include "regression/vli_test.hpp"

#include "vli/detail/bit_masks.hpp"
#include "tools.h"

#define Size1 3
#define Size2 2*Size1
#define Order 10

using vli::vli_cpu;
using vli::max_int_value;
using vli::monomial;
using vli::polynomial;
using vli::vector_polynomial;

typedef vli_cpu< unsigned long int, Size1> vli_type_cpu;
typedef vli_cpu< unsigned long int, Size2> vli_result_type_cpu;


//typedef vli::polynomial< vli_type_cpu, vli::max_order_each<Order>, vli::var<'x'> > polynomial_type_cpu;
//typedef vli::polynomial<vli_result_type_cpu  , vli::max_order_each<2*Order>, vli::var<'x'> > polynomial_result_type_cpu;

//typedef vli::polynomial< vli_type_cpu, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'> > polynomial_type_cpu;
//typedef vli::polynomial< vli_result_type_cpu, vli::max_order_each<2*Order>, vli::var<'x'>, vli::var<'y'> > polynomial_result_type_cpu;

typedef vli::polynomial< vli_type_cpu, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_type_combined_xyzw;
typedef vli::polynomial< vli_result_type_cpu, vli::max_order_combined<2*Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_result_type_combined_xyzw;

//typedef vli::polynomial< vli_type_cpu, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'> > polynomial_type_combined_cpu;
//typedef vli::polynomial< vli_result_type_cpu, vli::max_order_combined<2*Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'> > polynomial_result_type_combined_cpu;

typedef vli::polynomial< vli_type_cpu, vli::max_order_each<Order>, vli::var<'x'>  >polynomial_type_each_x;
typedef vli::polynomial< vli_result_type_cpu, vli::max_order_each<2*Order>, vli::var<'x'> > polynomial_result_type_each_x;

//typedef vli::polynomial< vli_type_cpu, vli::max_order_each<Order>, vli::var<'x'> , vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_type_cpu;
//typedef vli::polynomial< vli_result_type_cpu, vli::max_order_each<2*Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_result_type_cpu;

//typedef vli::vector_polynomial<polynomial_type_cpu> vector_type_cpu;
//typedef vli::vector_polynomial<polynomial_type_combined_cpu> vector_type_combined_cpu;

/*
int main (int argc, char * const argv[]) 
{

 
    int SizeVector = atoi(argv[1]);   

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
 polynomial_type_cpu tooo;
 polynomial_result_type_cpu result_pure_cpu  ;
 polynomial_result_type_cpu result_pure_cpu_omp  ;
 polynomial_result_type_cpu result_mix_cpu_gpu  ;
 polynomial_result_type_cpu result_cpu_gpu  ;
  

polynomial_type_combined_cpu polyc;
polynomial_result_type_combined_cpu polycres, polycres2;


 vector_type_combined_cpu v1c(SizeVector);
 vector_type_combined_cpu v2c(SizeVector);

    tools::fill_vector_random(v1c);
    tools::fill_vector_random(v2c);

    polycres2= vli::detail::inner_product_cpu(v1c,v2c);
    polycres = vli::detail::inner_product_gpu(v1c,v2c);

     if(polycres2 == polycres ) {printf("OK gpu combined\n"); } else{printf("NO OK gpu combined \n"); } 
    tools::fill_vector_random(v1);
    tools::fill_vector_random(v2);

    Timer t3("CPU vli_omp");
    t3.begin();
      result_pure_cpu = vli::detail::inner_product_cpu(v1,v2);
    t3.end();
#ifdef VLI_USE_GPU
std::cout << " --------------------------- " << std::endl;
    Timer t("GPU omp 1" );
    t.begin();   
      result_pure_cpu_omp = vli::inner_product(v1,v2);
    t.end();

std::cout << " --------------------------- " << std::endl;
    Timer t1("GPU omp 2 ");
    t1.begin();   
      result_pure_cpu_omp = vli::inner_product(v1,v2);
    t1.end();
std::cout << " --------------------------- " << std::endl;
    Timer t2("GPU omp 3 ");
    t2.begin();   
      result_pure_cpu_omp = vli::inner_product(v1,v2);
    t2.end();


#endif

    Timer t4("CPU gmp_omp");
    t4.begin();
    pgmpd = inner_product(v1gmp,v2gmp);
    t4.end();

     if(result_pure_cpu_omp ==result_pure_cpu ) {printf("OK gpu\n"); } else{printf("NO OK gpu \n"); } 


    if(ValidatePolyVLI_PolyGMP(result_pure_cpu,pgmpd))
     std::cout << "validation GMP OK " << std::endl;

    return 0;
}

*/

    typedef boost::mpl::vector<polynomial_type_each_x, polynomial_type_combined_xyzw> polynomial_list;
   
    struct test_case {
   
       template <typename Polynomial>
       void operator()(Polynomial const&) {
             typedef typename vli::polynomial_multiply_result_type<Polynomial>::type Polynomial_res;
             typedef vli::vector_polynomial<Polynomial> vector_polynomial;
             typedef vli::vector_polynomial<Polynomial_res> vector_polynomial_res;

             vector_polynomial v1(128),v2(128);
             Polynomial_res p1_res, p2_res;

             tools::fill_vector_random(v1);
             tools::fill_vector_random(v2);

             p1_res = vli::detail::inner_product_cpu(v1,v2);
             p2_res = vli::detail::inner_product_gpu(v1,v2);
           // if(p1_res == p2_res) {printf("OK gpu \n"); } else{printf("NO OK \n"); } 
           }
       };
   
   
   int main() {
   
       boost::mpl::for_each<polynomial_list>(test_case());
   
       return 0;
   
   }
