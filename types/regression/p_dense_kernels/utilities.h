
#include <boost/static_assert.hpp> 

namespace Random{
   struct random {
       random(){};
       double operator()(){return drand48();} 
       int IntRd(){return rand();}
   };
}

static Random::random Rd; //one generator is enough so static

template <int n, int m> // n # of workgroup, 
struct size {
   BOOST_STATIC_ASSERT(n>0);
   BOOST_STATIC_ASSERT(n*ambient::traits::value > m);
   typedef ambient::traits::value_type dbl; // To template later
   enum {valuex = n*ambient::traits::value+m};// n is the number or work group, m how we resize
   enum {valuey = n*ambient::traits::value-m};// n is the number or work group, m how we resize
   enum {null = 0};// n is the number or work group, m how we resize
};

typedef ambient::dim2 dim;
typedef maquis::types::dense_matrix<ambient::traits::value_type> sMatrix;
typedef maquis::types::p_dense_matrix<ambient::traits::value_type> pMatrix;
typedef maquis::types::diagonal_matrix<ambient::traits::value_type> sDiagMatrix;
typedef maquis::types::p_diagonal_matrix<ambient::traits::value_type> pDiagMatrix;
typedef boost::mpl::list<size<2,0>,size<2,3>,size<3,0>,size<3,3>,size<3,-3>,size<5,7>,size<5,-7> > test_types; // prime numbers ...

struct caveats {
    caveats() {
        ambient::init(); // init ambient
        srand48(1); //init random generator
        srand(1); //init random generator
    }

    ~caveats() {   
        ambient::finalize();
    }
};

BOOST_GLOBAL_FIXTURE( caveats );
