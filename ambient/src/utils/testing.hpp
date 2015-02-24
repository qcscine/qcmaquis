#ifndef TESTING_TOOLS
#define TESTING_TOOLS

#include "ambient/ambient.hpp"
#include "ambient/container/numeric/matrix.hpp"
#ifdef AMBIENT_ALPS
#include "alps/numeric/matrix.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "ambient/container/numeric/bindings/alps.hpp"
#endif
#include "ambient/container/numeric/bindings/types.hpp"

#define CATCH_CONFIG_MAIN
#include "utils/catch.hpp"
#include "utils/timings.hpp"

#define TOLERANCE (double)1/0xF4240
#define REQUIRE_CLOSE(a,b) REQUIRE(ambient::utils::require_close(a,b,TOLERANCE))

#define TEST_IB 256
#define TEST_M TEST_IB
#define TEST_N TEST_IB

namespace ambient { namespace utils {

    struct random {
        random(){};
        double operator()(){return drand48();} 
        int IntRd(){return rand();}
    } Rd;

    template<typename T>
    bool require_close(T left, T right, double tolerance){
        auto safe_div = [](double a, double b){
            if(a == 0 || (b > 1 && a < b*std::numeric_limits<double>::epsilon())) return (double)0;
            if(b < 1 && a > b*std::numeric_limits<double>::max()) return std::numeric_limits<double>::max();
            return a / b;
        };
        double d = std::abs(left - right);
        double l = safe_div(d, std::abs(left));
        double r = safe_div(d, std::abs(right));
        return std::max(l,r) <= tolerance;
    }

    template<size_t P>
    struct scenario {
        scenario(int argc = 0, char** argv = NULL) : argc(argc), argv(argv) 
        {
        }
        size_t num_threads(){
            return have_input(1) ? atoi(argv[1]) : 1;
        }
        size_t num_cols(){
            return have_input(2) ? atoi(argv[2]) : P;
        }

        size_t num_rows(){
            return have_input(3) ? atoi(argv[3]) : num_cols();
        }
        template<class FLOPS>
        void report(FLOPS flops, double time){
               std::cout << "-------------------------\n"
                         << " Time     " << time                                << "\n"
                         << " GFlops   " << flops(num_cols(), num_rows(), time) << "\n" 
                         << " Threads: " << num_threads()                       << "\n"
                         << " Matrix:  " << num_rows() << "x" << num_cols()     << "\n"
                         << "-------------------------\n";
        }
    private:
        bool have_input(size_t field){
            return ((argc > field) ? true : false);
        }
        int argc;
        char** argv;
    };

    struct measurement : scenario<AMBIENT_IB> {
        typedef ambient::async_timer timer;
    };
    struct gflops {
        static double gemm(size_t x, size_t y, double time){
            return 2*(double)x*(double)y*(double)y/(time*1.0e9);
        };
    };

} }


template<class T> using matrix = ambient::numeric::tiles<ambient::numeric::matrix<T> >;
template<class T> using diagonal = ambient::numeric::tiles<ambient::numeric::diagonal_matrix<T> >;

#ifdef AMBIENT_ALPS
template<class T> using matrix_ = alps::numeric::matrix<T>;
template<class T> using diagonal_ = alps::numeric::diagonal_matrix<T>;
#endif

using namespace ambient::numeric::bindings;
using ambient::utils::measurement;
using ambient::utils::gflops;
#endif
