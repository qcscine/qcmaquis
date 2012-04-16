#include <boost/static_assert.hpp> 
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/facilities/empty.hpp>

#define MAX_ITERATION 1
#define WorkGroup 1024
#define N 0xF4240 
#define ValueType double

#define TIM_PP_COMMA_IF(n) BOOST_PP_COMMA_IF( BOOST_PP_SUB(MAX_ITERATION,BOOST_PP_ADD(n,1))) 

namespace Random{
   struct random {
       random(){};
       double operator()(){return drand48();} 
       int IntRd(){return rand();}
   };
}

static Random::random Rd; //one generator is enough so static
static double value = (double)1/N; // for boost check close

template <int n, int m, typename T, int nthreads> // n # of workgroup, T double or std::complex<double> 
struct size {
   BOOST_STATIC_ASSERT(n>0);
   typedef T value_type; // To template later
   enum {ValueX = n};// n is the number or work group, m how we resize
   enum {ValueY = m};// n is the number or work group, m how we resize
   enum {ValueThread = nthreads}; // number of threads for blas
   enum {null = 0};// n is the number or work group, m how we resize
};

/*
* size<A,B,C>
* A BOOST_PP_ADD(n,1)*WorkGroup, size of the matrix (square matrix for the bench)
* B double, the type of the matrix
* C 2, the num of thread for the gemm classic test
*/

typedef boost::mpl::list< 
                     //  size<1024,1042,double,4>  // basic example if you do it by hand
                        #define DECLARATION_test_cases(z, n, unused) \
                            size< BOOST_PP_ADD(n,1)*WorkGroup, BOOST_PP_ADD(n,1)*WorkGroup ,double, 8> TIM_PP_COMMA_IF(n)  
                            BOOST_PP_REPEAT(MAX_ITERATION, DECLARATION_test_cases, ~)
                        #undef DECLARATION_test_cases

                        > test_types; // close >

struct caveats {
    caveats() {
        srand48(1); //init random generator
        srand(1); //init random generator
    }

    ~caveats() {
    }
};

void Boost_check_close_adapter(double a, double b){ 
    BOOST_CHECK_CLOSE(a, b, value); 
};

void Boost_check_close_adapter(std::complex<double> a, std::complex<double> b){ 
    BOOST_CHECK_CLOSE(a.real(), b.real(), value); 
    BOOST_CHECK_CLOSE(a.imag(), b.imag(), value); 
};

double GFlopsGemm(int n, int m, int k, double time){
    return 2*(double)n*(double)m*(double)k/(time*1.0e9); //get GigaFlop/s
};

BOOST_GLOBAL_FIXTURE( caveats );
