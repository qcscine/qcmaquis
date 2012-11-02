#include <boost/mpl/for_each.hpp>
#include <gmpxx.h>
#include <iomanip> 
#ifdef VLI_USE_GPU
#include "vli/detail/gpu/inner_product_gpu_accelerator.hpp"
#endif //VLI_USE_GPU

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial.hpp"
#include "vli/vli.hpp"
#include <time.h>
#include "utils/timings.h"
#include "utils/tools.h"
//#define Size_vec 65535// play with this 1024 - 16384
#define Size_vec 4096// play with this 1024 - 16384
#define Order 10// play 5 - 15, cautious outside memory, xyzw poly ( 10 is the real target)

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

typedef boost::mpl::vector<
                            polynomial_type_each_xyz_128,
                            polynomial_type_each_xy_128,
                            polynomial_type_each_x_128
//                            polynomial_type_each_xyzw_128// buffer can be too large cpu/gpu, be cautious
                          > polynomial_list_128_each;

typedef boost::mpl::vector<
                            polynomial_type_combined_xyzw_128,
                            polynomial_type_combined_xyz_128,
                            polynomial_type_combined_xy_128,
                            polynomial_type_combined_x_128 
                          > polynomial_list_128_combined;

typedef boost::mpl::vector<
//                         polynomial_type_each_xyzw_192
                           polynomial_type_each_xyz_192,
                           polynomial_type_each_xy_192,
                           polynomial_type_each_x_192 
                          > polynomial_list_192_each;

typedef boost::mpl::vector<
                           polynomial_type_combined_xyzw_192,
                           polynomial_type_combined_xyz_192,
                           polynomial_type_combined_xy_192,
                           polynomial_type_combined_x_192
                          > polynomial_list_192_combined;

typedef boost::mpl::vector<
                       //   polynomial_type_each_xyzw_256// buffer can be too large cpu/gpu, be cautious
                            polynomial_type_each_xyz_256,
                            polynomial_type_each_xy_256,
                            polynomial_type_each_x_256
                          > polynomial_list_256_each;

typedef boost::mpl::vector<
                            polynomial_type_combined_xyzw_256,
                            polynomial_type_combined_xyz_256,
                            polynomial_type_combined_xy_256,
                            polynomial_type_combined_x_256
                          > polynomial_list_256_combined;
   struct test_case {
    

   template <typename Polynomial>
   void operator()(Polynomial const&) {
       std::cout.precision(5);
       //GMP polys give by class traits
       typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type Polynomial_gmp;
       typedef typename vli::polynomial_multiply_type_gmp<Polynomial>::type_res Polynomial_gmp_res;
       typedef vli::vector_polynomial<Polynomial_gmp> vector_polynomial_gmp;
       typedef vli::vector_polynomial<Polynomial_gmp_res> vector_polynomial_gmp_res;
       //VLI poly
       typedef typename vli::polynomial_multiply_result_type<Polynomial>::type Polynomial_res;
       typedef vli::vector_polynomial<Polynomial> vector_polynomial;
       typedef vli::vector_polynomial<Polynomial_res> vector_polynomial_res;
       { 
       //VLI polys
       vector_polynomial v1(Size_vec),v2(Size_vec);
       Polynomial_res p1_res, p2_res;
       //GMP polys
       vector_polynomial_gmp v1_gmp(Size_vec),v2_gmp(Size_vec);
       Polynomial_gmp_res p_gmp_res;
        
       tools::fill_vector_random(v1);
       tools::fill_vector_random(v2);

       tools::converter(v1,v1_gmp); 
       tools::converter(v2,v2_gmp); 

       Timer tgmp("CPU GMP ");
       tgmp.begin();
           p_gmp_res = vli::detail::inner_product_cpu(v1_gmp,v2_gmp);
       tgmp.end();
        
       Timer t0("CPU ");
       t0.begin();
           p1_res = vli::detail::inner_product_cpu(v1,v2);
       t0.end();

       #ifdef VLI_USE_GPU
       Timer t1("GPU ");
       t1.begin();
           p2_res =  vli::detail::inner_product_gpu_helper<Polynomial>::inner_product_gpu(v1,v2);
       t1.end();

       if(tools::equal<Polynomial>(p1_res,p_gmp_res))
           if(p1_res == p2_res){ 
               std::cout << " ----- OK, cpu ";
               std::cout << t0.get_time() << " gpu ";
               std::cout << t1.get_time() << " gmp ";
               std::cout << tgmp.get_time();
               std::cout.precision(2);
               std::cout << " G vli: "  << tgmp.get_time()/t0.get_time() << " G gpu: " << tgmp.get_time()/t1.get_time()   ; 
           }else
               std::cout << " NO OK, " ; 
       #endif
       }
      

       }
   };

namespace vli { namespace detail { 
   void ultimate_192( boost::uint64_t* x/* %%rdi */,  boost::uint64_t const* y/* %%rsi */,  boost::uint64_t const* z/* %%rdx -> rbx */);
}}

int main(int argc, char* argv[]) {
       std::cout << " -------ASCII ART ^_^' --------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -------Size vector : " << Size_vec  << ", Order " << Order << std::endl;
       std::cout << " -----  MaxOrderEach------------------------------------------------------------------------------------------------------------------------------------------------------------ " << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  3 variable --------------------------------- 2 variables --------------------------------- 1 variables ----------------------------------------------------------------- " << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  256bits * 256bits = 512 bits ------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       boost::mpl::for_each<polynomial_list_256_each>(test_case());
       std::cout << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  192bits * 192bits = 384 bits ------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       boost::mpl::for_each<polynomial_list_192_each>(test_case());
       std::cout << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  128bits * 128bits = 256 bits ------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       boost::mpl::for_each<polynomial_list_128_each>(test_case());
       std::cout << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  MaxOrderCombined-------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  4 variable --------------------------------- 3 variables --------------------------------- 2 variables --------------------------------- 1 variables ------------------- " << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  256bits * 256bits = 512 bits ------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       boost::mpl::for_each<polynomial_list_256_combined>(test_case());
       std::cout << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  192bits * 192bits = 384 bits ------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       boost::mpl::for_each<polynomial_list_192_combined>(test_case());
       std::cout << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       std::cout << " -----  128bits * 128bits = 256 bits ------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       boost::mpl::for_each<polynomial_list_128_combined>(test_case());
       std::cout << std::endl;
       std::cout << " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- " << std::endl;
       return 0;
}
