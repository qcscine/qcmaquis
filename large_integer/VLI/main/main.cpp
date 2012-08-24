#include <boost/mpl/for_each.hpp>
 
#ifdef VLI_USE_GPU
#include "vli/detail/gpu/inner_product_gpu_booster.hpp"
//#include "vli/utils/gpu_manager.h"
//#include "vli/utils/gpu_manager.hpp"
#endif //VLI_USE_GPU

//#include "vli/polynomial/vector_polynomial_cpu.hpp"
//#include "vli/polynomial/polynomial.hpp"
#include "vli/vli.hpp"
#include "utils/timings.h"

//#include "tools.h"

#define Size_vec 1024
#define Order 10
#define Size_vli 3

using vlilib::vli;
//using vlilib::polynomial;
//using vlilib::vector_polynomial;
//typedef vli
typedef vli< 192> vli_type_cpu;
//typedef poly max order each
/*
typedef vli::polynomial< vli_type_cpu, vli::max_order_each<Order>, vli::var<'x'>  >polynomial_type_each_x;
typedef vli::polynomial< vli_type_cpu, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>  >polynomial_type_each_xy;
typedef vli::polynomial< vli_type_cpu, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>  >polynomial_type_each_xyz;
typedef vli::polynomial< vli_type_cpu, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'>  >polynomial_type_each_xyzw;
//typedef poly max order combined
typedef vli::polynomial< vli_type_cpu, vli::max_order_combined<Order>, vli::var<'x'> > polynomial_type_combined_x;
typedef vli::polynomial< vli_type_cpu, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'> > polynomial_type_combined_xy;
typedef vli::polynomial< vli_type_cpu, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'> > polynomial_type_combined_xyz;
typedef vli::polynomial< vli_type_cpu, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_type_combined_xyzw;

typedef boost::mpl::vector<polynomial_type_each_x,
                           polynomial_type_each_xy,
                           polynomial_type_each_xyz,
                         //  polynomial_type_each_xyzw  // buffer too large cpu/gpu
                           polynomial_type_combined_x,
                           polynomial_type_combined_xy,
                           polynomial_type_combined_xyz,
                           polynomial_type_combined_xyzw
                          > polynomial_list;

   struct test_case {

   template <typename Polynomial>
   void operator()(Polynomial const&) {
       typedef typename vli::polynomial_multiply_result_type<Polynomial>::type Polynomial_res;
       typedef vli::vector_polynomial<Polynomial> vector_polynomial;
       typedef vli::vector_polynomial<Polynomial_res> vector_polynomial_res;
       { 
       vector_polynomial v1(Size_vec),v2(Size_vec);
       Polynomial_res p1_res, p2_res;

       tools::fill_vector_random(v1);
       tools::fill_vector_random(v2);
      
       Timer t0("CPU omp");
       t0.begin();
       p1_res = vli::detail::inner_product_cpu(v1,v2);
       t0.end();
    
       #ifdef VLI_USE_GPU
       Timer t1("GPU ");
       t1.begin();
    //   p2_res =  vli::detail::inner_product_gpu_helper<Polynomial>::inner_product_gpu(v1,v2);
       t1.end();
      
       if(p1_res == p2_res) 
           printf("OK gpu \n"); 
       else
           printf("NO OK gpu \n");
       #endif
       }
       }
   };
   */
int main(int argc, char* argv[]) {
    std::cout << std::numeric_limits<unsigned long long int>::digits << std::endl;
    vli_type_cpu a(1);
    vli_type_cpu b(0);
    
    if( 0 >= b)
        std::cout << "ok" << std::endl;
    
    
    /*
    Timer A("tpotot");
    A.begin();
    
    for( long int i =0; i < 0xffffffff; ++i){
        vli_type_cpu a(-26+i);
        vli_type_cpu b(-2116+i);
        if( a != b){
        int c = 2;
            c += 1;}
    }
    
    A.end();
    */
//       boost::mpl::for_each<polynomial_list>(test_case());
       return 0;
   }
