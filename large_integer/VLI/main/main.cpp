#include <boost/mpl/for_each.hpp>
 
#ifdef VLI_USE_GPU
#include "vli/detail/gpu/inner_product_gpu_accelerator.hpp"
#endif //VLI_USE_GPU

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial.hpp"
#include "vli/vli.hpp"

#include "utils/timings.h"
#include "utils/tools.h"

#define Size_vec 128// play with this 1024 - 16384
#define Order 5 // play 5 - 15, cautious outside memory, xyzw poly ( 10 is the real target)

using vli::polynomial;
using vli::vector_polynomial;

typedef vli::vli<128> vli_type_cpu_128;
typedef vli::vli<192> vli_type_cpu_192;
typedef vli::vli<256> vli_type_cpu_256;
typedef vli::vli<384> vli_type_cpu_384;

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
// to check :  g++ -E -P -I /BOOST_PATH/include/ -I ../.. vli_number_cpu_function_hooks.hpp | sed  "s/n/;\\`echo -e '\n\r'`/g"  
namespace vli{
    namespace detail{
 void tutu( boost::uint64_t* x/* %%rdi */,  boost::uint64_t const* y/* %%rsi */,  boost::uint64_t const* z/* %%rdx -> rbx */);
 void bubu( boost::uint64_t* x/* %%rdi */,  boost::uint64_t const* y/* %%rsi */,  boost::uint64_t const* z/* %%rdx -> rbx */);
   }}
int main(int argc, char* argv[]) {
       vli_type_cpu_192 a,b;
//       vli_type_cpu_384 c,d;

      a[0] = 0xfff;
    //  a[1] = 0xfff;
    //  a[2] = 0xfff;
    
      vli_type_cpu_192 c(a), d(a);
    
      b[0] = 0xf;
   //   b[1] = 0xfff;
    std::cout   << a << std::endl;
    std::cout   << b << std::endl;

    
    vli_type_cpu_192 f = c % b;
    
      std::cout   << f << std::endl;
      c /= b;
      d %= b;
    
       std::cout << std::hex  << a << std::endl;
       std::cout << std::hex  << (d+(c*=b)) << std::endl;
       
    
//       a *= 4;
       
       
 //      std::cout <<  a << std::endl;
        
/*
       a[0] = 0xfffffff;
       a[1] = 0xfffffff;
       a[2] = 0;
       b[0] = 0xfffffff;
       b[1] = 0xfffffff    ;
       b[2] = 0;
       
       Timer t0("new");
       Timer t1("old");
     
       t1.begin();
       for(int i=0; i < 0xffffff; ++i)
           vli::detail::mul<3>(&d[0],&a[0],&b[0]);
       t1.end();

       t0.begin();
       for(int i=0; i < 0xffffff; ++i)
           vli::detail::bubu(&c[0],&a[0],&b[0]);
       t0.end();


       std::cout << std::hex << c << std::endl;
       std::cout << std::hex << d << std::endl;
     
       if(d == c)
           std::cout << " OK " << std::endl;

       boost::mpl::for_each<polynomial_list_128>(test_case());
       boost::mpl::for_each<polynomial_list_192>(test_case());
       boost::mpl::for_each<polynomial_list_256>(test_case());
       */
       return 0;
}
