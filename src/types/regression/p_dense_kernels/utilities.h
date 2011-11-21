
#include <boost/static_assert.hpp> 

namespace Random{
   struct random {
       random(){};
       double operator()(){return drand48();} 
       int IntRd(){return mrand48();}
   };
}

static Random::random Rd; //one generator is enough so static

template <int n, int m> // n # of workgroup, 
struct size {
   BOOST_STATIC_ASSERT(n>0);
   BOOST_STATIC_ASSERT(n*ambient::traits::value > m);
   typedef double dbl; // To template later
   enum {valuex = n*ambient::traits::value+m};// n is the number or work group, m how we resize
   enum {valuey = n*ambient::traits::value-m};// n is the number or work group, m how we resize
};

typedef ambient::dim2 dim;
typedef maquis::types::dense_matrix<double> sMatrix;
typedef maquis::types::p_dense_matrix<double> pMatrix;
typedef maquis::types::diagonal_matrix<double> sDiagMatrix;
typedef maquis::types::p_diagonal_matrix<double> pDiagMatrix;
typedef boost::mpl::list<size<2,0>,size<2,3>,size<3,0>,size<3,3>,size<3,-3>,size<5,7>, size<5,-7>  > test_types; // prime number ...

struct caveats {
    caveats() {
        ambient::init(); // init ambient
        srand48(1); //init random generator
    }

    ~caveats() {   
        ambient::finalize();
    }
};

BOOST_GLOBAL_FIXTURE( caveats );
