#include "vli/utils/macro.h"

namespace vli{
    namespace detail{
                     // new functions type : VLI<n*64> + VLI<n*64> : add128_128, add192_192 ...
                     #define FUNCTION_add_nbits_nbits(z, n, unused) \
                         void NAME_ADD_NBITS_PLUS_NBITS(n)(unsigned long int* x, unsigned long int const* y){ \
                         asm(                                                                                 \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), LOAD_register, ~)                         \
                                 "addq (%%rsi), "R(0)" \n"                                                    \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,1), ADC_register, ~)                          \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), SAVE_register, ~)                         \
                                 : : :BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), CLOTHER_register, ~) "memory"        \
                            );                                                                                \
                         }                                                                                    \

                     BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_add_nbits_nbits, ~)
                     #undef FUNCTION_add_nbits_nbits

                     //new functions type : VLI<n*64> + VLI<64> : add192_64, add256_64
                     //the case is done after add128_64
                     #define FUNCTION_add_nbits_64bits(z, n, unused) \
                         void NAME_ADD_NBITS_PLUS_64BITS(n)(unsigned long int* x, unsigned long int const* y){ \
                         asm(                                                                                  \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), LOAD_register, ~)                          \
                                 "addq (%%rsi), "R(0)" \n"                                                     \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,1), ADC0_register, ~)                          \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), SAVE_register, ~)                          \
                                 : : :BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), CLOTHER_register, ~) "memory"         \
                            );                                                                                 \
                         }                                                                                     \

                     BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_add_nbits_64bits, ~)
                     #undef FUNCTION_add_nbits_64bits

                     //new functions type : VLI<n*64> + VLI<(n-1)*64> : add128_64, add192_128 ...
                     #define FUNCTION_add_nbits_nminus1bits(z, n, unused) \
                         void NAME_ADD_NBITS_PLUS_NMINUS1BITS(n)(unsigned long int* x, unsigned long int const* y){ \
                         asm(                                                                                       \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), LOAD_register, ~)                               \
                                 "addq (%%rsi), "R(0)" \n"                                                          \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,1), ADC_register, ~)                                \
                                 "adcq $0x0, "R(BOOST_PP_ADD(n,2))" \n"                                             \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,3), SAVE_register, ~)                          \
                                 : : :BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), CLOTHER_register, ~) "memory"              \
                            );                                                                                      \
                        }                                                                                           \

                     BOOST_PP_REPEAT(MAX_ITERATION_MINUS_ONE, FUNCTION_add_nbits_nminus1bits, ~)
                     #undef FUNCTION_add_nbits_nminus1bits

                    } // end namespace detail
             } // end namespace vli
