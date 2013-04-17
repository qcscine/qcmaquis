#ifndef TESTING_TOOLS
#define TESTING_TOOLS

#include <mpi.h>
#include <omp.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/random.hpp>
#include <boost/static_assert.hpp> 

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utils/bindings.hpp"

#define BOOST_CLOSE (double)1/0xF4240

template <int n, int m, typename T, int t>
struct input{
   BOOST_STATIC_ASSERT(n > 0 && m > 0);
   typedef T value_type;
   enum {valuex = n};
   enum {valuey = m};
   enum {valuet = t};
};

struct caveats{
    caveats(){ srand48(1); srand(1); }
    ~caveats(){}
};

BOOST_GLOBAL_FIXTURE( caveats );

struct random {
    random(){};
    double operator()(){return drand48();} 
    int IntRd(){return rand();}
} Rd;


void Boost_check_close_adapter(double a, double b){ 
    BOOST_CHECK_CLOSE(a, b, BOOST_CLOSE); 
};

void Boost_check_close_adapter(std::complex<double> a, std::complex<double> b){ 
    BOOST_CHECK_CLOSE(a.real(), b.real(), BOOST_CLOSE); 
    BOOST_CHECK_CLOSE(a.imag(), b.imag(), BOOST_CLOSE); 
};


bool have_input(size_t field){
    return ((boost::unit_test::framework::master_test_suite().argc > field) ? true : false);
}

template<typename T>
size_t get_input_threads(){
    size_t threads = T::valuet;
    if(have_input(1)) threads = atoi(boost::unit_test::framework::master_test_suite().argv[1]);
    return threads;
}

template<typename T>
size_t get_input_x(){
    size_t dim = T::valuex;
    if(have_input(2)) dim = atoi(boost::unit_test::framework::master_test_suite().argv[2]);
    return dim;
}

template<typename T>
size_t get_input_y(){
    size_t dim = T::valuey;
    if(have_input(3)) dim = atoi(boost::unit_test::framework::master_test_suite().argv[3]);
    else if(have_input(2)) dim = atoi(boost::unit_test::framework::master_test_suite().argv[2]);
    return dim;
}

template<class TIMER>
void report(const TIMER& a, double(*gflops)(size_t, size_t, double), size_t x, size_t y, size_t nthreads){
       std::cout << "-------------------------\n"
                 << " Time     " << a.get_time()               << "\n"
                 << " GFlops   " << gflops(x,y,a.get_time())   << "\n" 
                 << " Threads: " << nthreads                   << "\n"
                 << " Matrix:  " << y << "x" << x              << "\n"
                 << "-------------------------\n";
}

double GFlopsGemm(size_t x, size_t y, double time){
    return 2*(double)x*(double)y*(double)y/(time*1.0e9);
};

using namespace maquis::bindings;
#endif
