#include <boost/mpl/for_each.hpp>
//#include <gmp.h>
 
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
#define Size_vec 16384// play with this 1024 - 16384
#define Order 10 // play 5 - 15, cautious outside memory, xyzw poly ( 10 is the real target)

using vli::polynomial;
using vli::vector_polynomial;
//vli::polynomial<mpq_class, max_order_each<Order>, var<'x'> >
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
//                            polynomial_type_each_xyz_128,
                 //           polynomial_type_each_xyzw_128,// buffer can be too large cpu/gpu, be cautious
                            polynomial_type_combined_x_128,
                            polynomial_type_combined_xy_128
 //                           polynomial_type_combined_xyz_128
                         //   polynomial_type_combined_xyzw_128// buffer can be too large cpu/gpu, be cautious
                          > polynomial_list_128;

typedef boost::mpl::vector<// polynomial_type_each_x_192,
                            polynomial_type_each_xy_192
  //                          polynomial_type_each_xyz_192,
                          //  polynomial_type_each_xyzw_192,// buffer can be too large cpu/gpu, be cautious
                        //    polynomial_type_combined_x_192,
                       //     polynomial_type_combined_xy_192
//                            polynomial_type_combined_xyz_192
                         //   polynomial_type_combined_xyzw_192// buffer can be too large cpu/gpu, be cautious
                          > polynomial_list_192;

typedef boost::mpl::vector< polynomial_type_each_x_256,
                            polynomial_type_each_xy_256,
                            polynomial_type_each_xyz_256,
             //               polynomial_type_each_xyzw_256,// buffer can be too large cpu/gpu, be cautious 
                            polynomial_type_combined_x_256,
                            polynomial_type_combined_xy_256,
                            polynomial_type_combined_xyz_256
//                            polynomial_type_combined_xyzw_256// buffer can be too large cpu/gpu, be cautious
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
           std::cout << "OK gpu " << std::endl; 
       else
           std::cout << "NO OK gpu " << std::endl; 
    
       #endif
       }
      

       }
   };

namespace vli { namespace detail { 
   void ultimate_192( boost::uint64_t* x/* %%rdi */,  boost::uint64_t const* y/* %%rsi */,  boost::uint64_t const* z/* %%rdx -> rbx */);
}}

int main(int argc, char* argv[]) {
       vli_type_cpu_192 a,b;
       vli_type_cpu_384 c,d;   

       a[0] = 1;
       b[0] = 1;

//       a[0]=325165457354384;
//       a[1]=135354587344384;
//       a[2]=554984684684688;
//
//       b[0]=325165457354384;
//       b[1]=135354587344384;
//       b[2]=554984684684688;

       std::cout << "a :" << std::hex << a <<std::endl;
       std::cout << "b :" << std::hex << b <<std::endl;

       std::cout << std::endl << "+,+" <<std::endl;
       vli::detail::ultimate_192(&c[0],&a[0],&b[0]); 
       std::cout <<"ultimate: " << std::hex << c << std::endl;
       vli::detail::mul<3>(&d[0],&a[0],&b[0]);
       std::cout <<"correct:  " << std::hex << d << std::endl;

       b.negate();
       std::cout << std::endl << "+,-" <<std::endl;
       vli::detail::ultimate_192(&c[0],&a[0],&b[0]); 
       std::cout <<"ultimate: " << std::hex << c << std::endl;
       vli::detail::mul<3>(&d[0],&a[0],&b[0]);
       std::cout <<"correct:  " << std::hex << d << std::endl;

       a.negate();
       b.negate();
       std::cout << std::endl << "-,+" <<std::endl;
       vli::detail::ultimate_192(&c[0],&a[0],&b[0]); 
       std::cout <<"ultimate: " << std::hex << c << std::endl;
       vli::detail::mul<3>(&d[0],&a[0],&b[0]);
       std::cout <<"correct:  " << std::hex << d << std::endl;

       b.negate();
       std::cout << std::endl << "-,-" <<std::endl;
       vli::detail::ultimate_192(&c[0],&a[0],&b[0]); 
       std::cout <<"ultimate: " << std::hex << c << std::endl;
       vli::detail::mul<3>(&d[0],&a[0],&b[0]);
       std::cout <<"correct:  " << std::hex << d << std::endl;


       timespec ta,tb,tc;
       clock_gettime(CLOCK_MONOTONIC,&ta);
       for(unsigned int i=0; i < 100000000; ++i)
           vli::detail::ultimate_192(&c[0],&a[0],&b[0]); 

       clock_gettime(CLOCK_MONOTONIC,&tb);
       for(unsigned int i=0; i < 100000000; ++i)
           vli::detail::mul<3>(&c[0],&a[0],&b[0]); 
       clock_gettime(CLOCK_MONOTONIC,&tc);


       std::cout << "ultimate t=" << (tb.tv_sec - ta.tv_sec) + 1e-9 * (tb.tv_nsec - ta.tv_nsec) << "s" << std::endl;
       std::cout << "old mul  t=" << (tc.tv_sec - tb.tv_sec) + 1e-9 * (tc.tv_nsec - tb.tv_nsec) << "s" << std::endl;
       clock_gettime(CLOCK_MONOTONIC,&ta);
       for(unsigned int i=0; i < 100000000; ++i)
           vli::detail::mul<3>(&c[0],&a[0],&b[0]); 
       clock_gettime(CLOCK_MONOTONIC,&tb);
       for(unsigned int i=0; i < 100000000; ++i)
           vli::detail::ultimate_192(&c[0],&a[0],&b[0]); 
       clock_gettime(CLOCK_MONOTONIC,&tc);


       std::cout << "old mul  t=" << (tb.tv_sec - ta.tv_sec) + 1e-9 * (tb.tv_nsec - ta.tv_nsec) << "s" << std::endl;
       std::cout << "ultimate t=" << (tc.tv_sec - tb.tv_sec) + 1e-9 * (tc.tv_nsec - tb.tv_nsec) << "s" << std::endl;

//       boost::mpl::for_each<polynomial_list_128>(test_case());
//       boost::mpl::for_each<polynomial_list_192>(test_case());
//       boost::mpl::for_each<polynomial_list_256>(test_case());

       return 0;
}
