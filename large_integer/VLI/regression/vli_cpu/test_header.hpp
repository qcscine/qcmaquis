#ifndef VLI_TEST_HEADER_HPP
#define VLI_TEST_HEADER_HPP

#define BOOST_TEST_MODULE vli_cpu
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/preprocessor/cat.hpp>
#include <stdexcept>
#include <vli/vli_cpu.h>
#ifdef VLI_FUZZ_TESTS
#include <ctime>
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#endif // VLI_FUZZ_TESTS

#ifdef VLI_FUZZ_TESTS
#define VLI_FUZZABLE_TEST( TEST_NAME )                          \
    template <typename Initalizer>                              \
    static void BOOST_PP_CAT(TEST_NAME,_test(Initalizer init)); \
    BOOST_AUTO_TEST_CASE( TEST_NAME )                           \
    {                                                           \
        std::cout<<"Running "<<vli::test::fuzz_initializer::fuzz_iterations<<" iterations..."<<std::endl;    \
        for(unsigned int i=0; i < vli::test::fuzz_initializer::fuzz_iterations; ++i)                                  \
            BOOST_PP_CAT(TEST_NAME,_test(fuzz_initializer())); \
    }                                                           \
    template <typename Initalizer>                              \
    static void BOOST_PP_CAT(TEST_NAME,_test(Initalizer init))
#define VLI_STATIC_TEST( TEST_NAME )                            \
    template <typename Initalizer>                              \
    static void BOOST_PP_CAT(TEST_NAME,_test(Initalizer init))
#else  // VLI_FUZZ_TESTS
#define VLI_FUZZABLE_TEST( TEST_NAME ) VLI_STATIC_TEST( TEST_NAME )
#define VLI_STATIC_TEST( TEST_NAME )                            \
    template <typename Initalizer>                              \
    static void BOOST_PP_CAT(TEST_NAME,_test(Initalizer init)); \
    BOOST_AUTO_TEST_CASE( TEST_NAME )                           \
    {                                                           \
        BOOST_PP_CAT(TEST_NAME,_test(initializer()));           \
    }                                                           \
    template <typename Initalizer>                              \
    static void BOOST_PP_CAT(TEST_NAME,_test(Initalizer init))
#endif // VLI_FUZZ_TESTS

namespace vli {
namespace test {

typedef vli::vli_cpu<unsigned long int, VLI_SIZE> vli_type;
typedef vli::vli_cpu<unsigned long int, 2*VLI_SIZE> vli_type_double;

    
enum variant_enum {
      max_positive = 0      // fill with the max positive number
    , overflow_safe = 1     // fill such that x+x doesn't cause an overflow
    , fill_ff = 2           // fill with 0xff..ff
    , multiplies_overflow_safe = 3
};

struct initializer {
    void operator()(vli_type& v, variant_enum variant = overflow_safe) {
            for(std::size_t i=0; i != vli_type::size; ++i)
                v[i] = std::numeric_limits<vli_type::value_type>::max();
        switch(variant) {
            case overflow_safe:
                v[vli_type::size-1] = 1;
                break;
            case max_positive:
                v[vli_type::size-1] = 0x7fffffffffffffff;
                break;
            case multiplies_overflow_safe:
                v[vli_type::size/2-1] = 0x7fffffffffffffff;
                for(std::size_t i=vli_type::size/2; i != vli_type::size; ++i)
                    v[i] = 0;
                break;
            case fill_ff:
                break;
            default:
                throw(std::runtime_error("init variant not implemented"));
        }
    }
    void operator()(vli_type_double& v, variant_enum variant = overflow_safe) {
        for(std::size_t i=0; i != vli_type_double::size; ++i)
            v[i] = std::numeric_limits<vli_type::value_type>::max();
        switch(variant) {
            case max_positive:
                v[vli_type_double::size-1] = v[vli_type_double::size-1] & ~(vli_type::value_type(1)<<(sizeof(vli_type_double::value_type)*8-1));
                break;
            case overflow_safe:
                v[vli_type_double::size-1] = 1;
                break;
            case multiplies_overflow_safe:
                v[vli_type_double::size/2-1] = v[vli_type_double::size/2-1] & 0x7fffffffffffffff;
                for(std::size_t i=vli_type_double::size/2; i != vli_type_double::size; ++i)
                    v[i] = 0;
                break;
            default:
                throw(std::runtime_error("init variant not implemented"));
       }
    }
    void operator()(int& i, variant_enum variant = overflow_safe) {
        switch(variant) {
            case overflow_safe:
                i = std::numeric_limits<int>::max()/2;
                break;
            case max_positive:
                i = std::numeric_limits<int>::max();
                break;
            case fill_ff:
                i = -1;
                break;
            default:
                throw(std::runtime_error("init variant not implemented"));
        }
    }
};
#ifdef VLI_FUZZ_TESTS
struct fuzz_initializer {
      static boost::random::mt19937                                           rng;
      static boost::random::uniform_int_distribution<vli_type::value_type>    vli_value_type_max_rnd;
      static boost::random::uniform_int_distribution<int>                     int_plus_rnd;
    static unsigned int                                                     fuzz_iterations;
    void operator()(vli_type& v, variant_enum variant = overflow_safe) {
        for(std::size_t i=0; i != vli_type::size; ++i)
            v[i] = vli_value_type_max_rnd(rng);
        switch(variant) {
            case max_positive:
                v[vli_type::size-1] = v[vli_type::size-1] & ~(vli_type::value_type(1)<<(sizeof(vli_type::value_type)*8-1));
                break;
            case overflow_safe:
                v[vli_type::size-1] = 1;
                break;
            case multiplies_overflow_safe:
                v[vli_type::size/2-1] = v[vli_type::size/2-1] & 0x7fffffffffffffff;
                for(std::size_t i=vli_type::size/2; i != vli_type::size; ++i)
                    v[i] = 0;
                break;
            default:
                throw(std::runtime_error("init variant not implemented"));
        }
    }
    void operator()(vli_type_double& v, variant_enum variant = overflow_safe) {
        for(std::size_t i=0; i != vli_type_double::size; ++i)
            v[i] = vli_value_type_max_rnd(rng);
        switch(variant) {
            case max_positive:
                v[vli_type_double::size-1] = v[vli_type_double::size-1] & ~(vli_type::value_type(1)<<(sizeof(vli_type_double::value_type)*8-1));
                break;
            case overflow_safe:
                v[vli_type_double::size-1] = 1;
                break;
            case multiplies_overflow_safe:
                v[vli_type_double::size/2-1] = v[vli_type_double::size/2-1] & 0x7fffffffffffffff;
                for(std::size_t i=vli_type_double::size/2; i != vli_type_double::size; ++i)
                    v[i] = 0;
                break;
            default:
                throw(std::runtime_error("init variant not implemented"));
        }
    }
    void operator()(int& i, variant_enum variant = overflow_safe) {
        i = int_plus_rnd(rng);
        switch(variant) {
            // TODO this will always give a positive number
            case overflow_safe:
                i/=2;
                break;
            default:
                throw(std::runtime_error("init variant not implemented"));
        }
    }
};
unsigned int                                                  fuzz_initializer::fuzz_iterations        = 0;
boost::random::mt19937                                        fuzz_initializer::rng;
boost::random::uniform_int_distribution<vli_type::value_type> fuzz_initializer::vli_value_type_max_rnd = boost::random::uniform_int_distribution<vli_type::value_type>(0,std::numeric_limits<vli_type::value_type>::max());
boost::random::uniform_int_distribution<int>                  fuzz_initializer::int_plus_rnd           = boost::random::uniform_int_distribution<int>(0,std::numeric_limits<int>::max());
#endif // VLI_FUZZ_TESTS


template <typename T>
struct extended {
};

template <typename T, std::size_t Size>
struct extended<vli::vli_cpu<T,Size> > {
    typedef vli::vli_cpu<T,Size+1> type;
};

template <typename T>
struct double_sized_vli {
};

template <typename T, std::size_t Size>
struct double_sized_vli<vli::vli_cpu<T,Size> > {
    typedef vli::vli_cpu<T,2*Size> type;
};

} // end namespace test
} // end namespace vli

int main(int argc, char* argv[])
{
#ifdef VLI_FUZZ_TESTS

    if(argc > 1 && argc <= 3) {
        vli::test::fuzz_initializer::fuzz_iterations = boost::lexical_cast<unsigned int>(argv[1]);
        unsigned int seed = std::time(0);
        if (argc == 3)
            seed = boost::lexical_cast<unsigned int>(argv[2]);
        std::cout<<"Using random seed: "<<seed<<std::endl;
       // vli::test::fuzz_initializer::rng.seed(seed);

    } else {
        std::cerr<<"Usage:"<<std::endl;
        std::cerr<<argv[0]<<" <iterations> [<random seed>]"<<std::endl;
        return -1;
    }
#endif // VLI_FUZZ_TESTS
    return ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
}

#endif // VLI_TEST_HEADER_HPP
