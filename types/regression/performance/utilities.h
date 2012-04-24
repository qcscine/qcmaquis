#include <boost/static_assert.hpp> 

struct caveats{
    caveats(){ srand48(1); srand(1); }
    ~caveats(){}
};

template <int n, int m, typename T, int nthreads>
struct input{
   BOOST_STATIC_ASSERT(n > 0);
   typedef T value_type;
   enum {ValueX = n};
   enum {ValueY = m};
   enum {ValueThread = nthreads};
};

double GFlopsGemm(size_t x, size_t y, double time){
    return 2*(double)x*(double)y*(double)y/(time*1.0e9);
};

bool have_input(size_t field){
    return ((boost::unit_test::framework::master_test_suite().argc > field) ? true : false);
}

template<typename T>
size_t get_input_threads(){
    size_t threads = T::ValueThread;
    if(have_input(1)) threads = atoi(boost::unit_test::framework::master_test_suite().argv[1]);
    return threads;
}

template<typename T>
size_t get_input_x(){
    size_t dim = T::ValueX;
    if(have_input(2)) dim = atoi(boost::unit_test::framework::master_test_suite().argv[2]);
    return dim;
}

template<typename T>
size_t get_input_y(){
    size_t dim = T::ValueY;
    if(have_input(3)) dim = atoi(boost::unit_test::framework::master_test_suite().argv[3]);
    else if(have_input(2)) dim = atoi(boost::unit_test::framework::master_test_suite().argv[2]);
    return dim;
}

void report(const Timer& a, double(*gflops)(size_t, size_t, double), size_t x, size_t y, size_t nthreads){
    maquis::cout << "-------------------------\n"
                 << " Time     " << a.get_time()               << "\n"
                 << " GFlops   " << gflops(x,y,a.get_time())   << "\n" 
                 << " Threads: " << nthreads                   << "\n"
                 << " Matrix:  " << y << "x" << x              << "\n"
                 << "-------------------------\n";
}

typedef boost::mpl::list< input<1024,1024,double,1> > test_types;
BOOST_GLOBAL_FIXTURE(caveats);
