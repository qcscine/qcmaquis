/* \cond I do not need this part in the doc*/
//system and boost
#include <boost/mpl/for_each.hpp>
#include "boost/lexical_cast.hpp"
#include <gmpxx.h>
#include <iomanip> 
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

//#include "main/inline_add.h"

#ifdef VLI_USE_GPU
#include "vli/detail/gpu/inner_product_gpu_accelerator.hpp"
#endif //VLI_USE_GPU
//vli
#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial.hpp"
#include "vli/vli.hpp"
#include "vli/detail/kernels_cpu.h"
//utils
#include "utils/timings.h"
#include "utils/tools.h"
#include "misc.hpp"

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/stringize.hpp>

typedef vli::integer<128> integer_type_cpu_128;
typedef vli::integer<192> integer_type_cpu_192;
typedef vli::integer<256> integer_type_cpu_256;

typedef boost::mpl::vector<
                            integer_type_cpu_128,
                            integer_type_cpu_192,
                            integer_type_cpu_256
                          > vli_list;

   template<int NumBits>
   struct timescheduler{
       static void save(double tcpu, double tgmp){
           std::string name("time_bench");
           name += boost::lexical_cast<std::string>(NumBits);
           std::ofstream os(name.c_str(),std::ios::app);
           os << " " << NumBits << " " << tcpu << " " << tgmp << std::endl;
           os.close();
       };
   };   

   struct test_mul {
   template <typename integer>
       void operator()(integer const&) {
           typedef vli::integer<2*integer::numbits> integer_res;
           
           integer a,b,c;
           integer_res c_res;
       
           tools::fill_random(a);
           tools::fill_random(b);
       
            mpz_class agmp(a), bgmp(b), cgmp;
           
            Timer tvli("vli");
            Timer tgmp("gmp");
        
            long long limit=1;

            while( limit != 1000000000){
                     tvli.begin();
                     for (int i=1; i <= limit; i++)
                         multiply_extend(c_res,a,b);
                     tvli.end();
                     
                     tgmp.begin();
                     for (int i=1; i <= limit; i++)
                         cgmp = agmp * bgmp;

                     tgmp.end();
                     std::cout << " mul " << integer::numbits << " " << limit << " "  << tvli.get_time() << " " << tgmp.get_time() << std::endl;
                     limit *= 10;
            }

           mpz_class res(c_res);
           
           if(res == cgmp){
             std::cout << "ok " << std::endl;
           }
       }
   };

   struct test_add {
   template <typename integer>
       void operator()(integer const&) {
           typedef vli::integer<2*integer::numbits> integer_res;
           
           integer a,b,c;
           integer_res c_res;
       
           tools::fill_random(a);
           tools::fill_random(b);
       
            mpz_class agmp(a), bgmp(b), cgmp;
           
            Timer tvli("vli");
            Timer tgmp("gmp");
        
            int limit=1;

            while( limit != 1000000000){
                     tvli.begin();
                     for (int i=1; i <= limit; i++)
                         c = a + b;
                     tvli.end();
                     
                     tgmp.begin();
                     for (int i=1; i <= limit; i++)
                         cgmp = agmp + bgmp;

                     tgmp.end();
                     std::cout << " add " << integer::numbits << " " << limit << " "  << tvli.get_time() << " " << tgmp.get_time() << std::endl;
                     limit *= 10;
            }

           mpz_class res(c);
           
           if(res == cgmp){
             std::cout << "ok " << std::endl;
           }
       }
   };

int main(int argc, char* argv[]) {
       boost::mpl::for_each<vli_list>(test_add());
       boost::mpl::for_each<vli_list>(test_mul());
       return 0;
}
