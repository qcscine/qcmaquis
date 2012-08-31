#include <boost/mpl/for_each.hpp>
 
#ifdef VLI_USE_GPU
#include "vli/detail/gpu/inner_product_gpu_accelerator.hpp"
#endif //VLI_USE_GPU

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial.hpp"
#include "vli/vli.hpp"

#include "utils/timings.h"
#include "utils/tools.h"

#define Size_vec 1024 // play with this 1024 - 16384
#define Order 5 // play 5 - 15, cautious outside memory, xyzw poly ( 10 is the real target)

using vli::polynomial;
using vli::vector_polynomial;

typedef vli::vli<128> vli_type_cpu_128;
typedef vli::vli<192> vli_type_cpu_192;
typedef vli::vli<256> vli_type_cpu_256;

/*  -------------------------------------------------------------------- 128 bits ---------------------------------------------------------------------------------- */
typedef vli::polynomial< vli_type_cpu_128, vli::max_order_each<Order>, vli::var<'x'>  >polynomial_type_each_x_128;
typedef vli::polynomial< vli_type_cpu_128, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>  >polynomial_type_each_xy_128;
typedef vli::polynomial< vli_type_cpu_128, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>  >polynomial_type_each_xyz_128;
typedef vli::polynomial< vli_type_cpu_128, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'>  >polynomial_type_each_xyzw_128;
typedef vli::polynomial< vli_type_cpu_128, vli::max_order_combined<Order>, vli::var<'x'> > polynomial_type_combined_x_128;
typedef vli::polynomial< vli_type_cpu_128, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'> > polynomial_type_combined_xy_128;
typedef vli::polynomial< vli_type_cpu_128, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'> > polynomial_type_combined_xyz_128;
typedef vli::polynomial< vli_type_cpu_128, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_type_combined_xyzw_128;
/*  -------------------------------------------------------------------- 192 bits ---------------------------------------------------------------------------------- */
typedef vli::polynomial< vli_type_cpu_192, vli::max_order_each<Order>, vli::var<'x'>  >polynomial_type_each_x_192;
typedef vli::polynomial< vli_type_cpu_192, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>  >polynomial_type_each_xy_192;
typedef vli::polynomial< vli_type_cpu_192, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>  >polynomial_type_each_xyz_192;
typedef vli::polynomial< vli_type_cpu_192, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'>  >polynomial_type_each_xyzw_192;
typedef vli::polynomial< vli_type_cpu_192, vli::max_order_combined<Order>, vli::var<'x'> > polynomial_type_combined_x_192;
typedef vli::polynomial< vli_type_cpu_192, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'> > polynomial_type_combined_xy_192;
typedef vli::polynomial< vli_type_cpu_192, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'> > polynomial_type_combined_xyz_192;
typedef vli::polynomial< vli_type_cpu_192, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_type_combined_xyzw_192;
/*  -------------------------------------------------------------------- 256 bits ---------------------------------------------------------------------------------- */
typedef vli::polynomial< vli_type_cpu_256, vli::max_order_each<Order>, vli::var<'x'>  >polynomial_type_each_x_256;
typedef vli::polynomial< vli_type_cpu_256, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>  >polynomial_type_each_xy_256;
typedef vli::polynomial< vli_type_cpu_256, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>  >polynomial_type_each_xyz_256;
typedef vli::polynomial< vli_type_cpu_256, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'>  >polynomial_type_each_xyzw_256;
typedef vli::polynomial< vli_type_cpu_256, vli::max_order_combined<Order>, vli::var<'x'> > polynomial_type_combined_x_256;
typedef vli::polynomial< vli_type_cpu_256, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'> > polynomial_type_combined_xy_256;
typedef vli::polynomial< vli_type_cpu_256, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'> > polynomial_type_combined_xyz_256;
typedef vli::polynomial< vli_type_cpu_256, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_type_combined_xyzw_256;

typedef boost::mpl::vector< polynomial_type_each_x_128,
                            polynomial_type_each_xy_128,
                            polynomial_type_each_xyz_128,
                            polynomial_type_each_xyzw_128,// buffer can be too large cpu/gpu, be cautious
                            polynomial_type_combined_x_128,
                            polynomial_type_combined_xy_128,
                            polynomial_type_combined_xyz_128,
                            polynomial_type_combined_xyzw_128// buffer can be too large cpu/gpu, be cautious
                          > polynomial_list_128;

typedef boost::mpl::vector< polynomial_type_each_x_192,
                            polynomial_type_each_xy_192,
                            polynomial_type_each_xyz_192,
                            polynomial_type_each_xyzw_192,// buffer can be too large cpu/gpu, be cautious
                            polynomial_type_combined_x_192,
                            polynomial_type_combined_xy_192,
                            polynomial_type_combined_xyz_192,
                            polynomial_type_combined_xyzw_192// buffer can be too large cpu/gpu, be cautious
                          > polynomial_list_192;

typedef boost::mpl::vector< polynomial_type_each_x_256,
                            polynomial_type_each_xy_256,
                            polynomial_type_each_xyz_256,
                            polynomial_type_each_xyzw_256,// buffer can be too large cpu/gpu, be cautious 
                            polynomial_type_combined_x_256,
                            polynomial_type_combined_xy_256,
                            polynomial_type_combined_xyz_256,
                            polynomial_type_combined_xyzw_256// buffer can be too large cpu/gpu, be cautious
                          > polynomial_list_256;

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
       p2_res =  vli::detail::inner_product_gpu_helper<Polynomial>::inner_product_gpu(v1,v2);
       t1.end();
      
       if(p1_res == p2_res) 
           printf("OK gpu \n"); 
       else
           printf("NO OK gpu \n");
       #endif
       }

       }
   };
   
int main(int argc, char* argv[]) {
       boost::mpl::for_each<polynomial_list_128>(test_case());
       boost::mpl::for_each<polynomial_list_192>(test_case());
       boost::mpl::for_each<polynomial_list_256>(test_case());
       return 0;
}
