#include "utils/bindings.hpp"
#include <boost/static_assert.hpp> 

#define M 0xF4240 
#define ValueWG  128

namespace Random{
   struct random {
       random(){};
       double operator()(){return drand48();} 
       int IntRd(){return rand();}
   };
}

static Random::random Rd; //one generator is enough so static
static double value = (double)1/M; // for boost check close

template <int n, int m, typename T> // n # of workgroup, T double or std::complex<double> 
struct size {
   BOOST_STATIC_ASSERT(n>0);
   BOOST_STATIC_ASSERT(n*ValueWG > m);
   typedef T value_type; 
   enum {valuex = n*ValueWG+m};// n is the number or work group, m how we resize
   enum {valuey = n*ValueWG-m};// n is the number or work group, m how we resize
   enum {null = 0};// n is the number or work group, m how we resize
};

typedef ambient::dim2 dim;

typedef boost::mpl::list< 
/*
                         size<2,0,double>,
                         size<2,0,std::complex<double> >,
                         size<2,3,double>,
                         size<2,3,std::complex<double> >,
                         size<2,-3,double>,
                         size<2,-3,std::complex<double> >,
                         size<5,0,double>,
                         size<5,0,std::complex<double> >,
                         size<5,-7,double>,
                         size<5,-7,std::complex<double> >,
*/
                         size<5,-7,double>,
                         size<5,-7,std::complex<double> >
                        > test_types; 

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

BOOST_GLOBAL_FIXTURE( caveats );

