#include "vli/utils/macro.h"
#include <cassert>
namespace vli{
    namespace detail{
                     // new functions type : VLI<n*64> - VLI<n*64> : sub128_128, sub192_192 ...
                     #define FUNCTION_sub_nbits_nbits(z, n, unused) \
                         void NAME_SUB_NBITS_MINUS_NBITS(n)(unsigned long int* x, unsigned long int const* y){ \
                         asm(                                                                                 \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), LOAD_register, ~)                         \
                                 "subq (%%rsi), "R(0)" \n"                                                    \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,1), SBB_register, ~)                          \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), SAVE_register, ~)                         \
                                 : : :BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), CLOTHER_register, ~) "memory"        \
                            );                                                                                \
                         }                                                                                    \

                     BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_sub_nbits_nbits, ~)
                     #undef FUNCTION_sub_nbits_nbits

                     //new functions type : VLI<n*64> - VLI<64> : sub192_64, sub256_64
                     //the case is done after sub128_64
                     #define FUNCTION_sub_nbits_64bits(z, n, unused) \
                         void NAME_SUB_NBITS_MINUS_64BITS(n)(unsigned long int* x, unsigned long int const* y){ \
                         asm(                                                                                  \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), LOAD_register, ~)                          \
                                 "subq (%%rsi), "R(0)" \n"                                                     \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,1), SBB0_register, ~)                          \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), SAVE_register, ~)                          \
                                 : : :BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), CLOTHER_register, ~) "memory"         \
                            );                                                                                 \
                         }                                                                                     \

                     BOOST_PP_REPEAT(MAX_ITERATION, FUNCTION_sub_nbits_64bits, ~)
                     #undef FUNCTION_sub_nbits_64bits

                     //new functions type : VLI<n*64> - VLI<(n-1)*64> : sub128_64, sub192_128 ...
                     #define FUNCTION_sub_nbits_nminus1bits(z, n, unused) \
                         void NAME_SUB_NBITS_MINUS_NMINUS1BITS(n)(unsigned long int* x, unsigned long int const* y){ \
                             assert(false); \
                         asm(                                                                                       \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), LOAD_register, ~)                               \
                                 "subq (%%rsi), "R(0)" \n"                                                          \
                                 BOOST_PP_REPEAT(BOOST_PP_ADD(n,1), SBB_register, ~)                                \
                                 "sbbq $0x0, "R(BOOST_PP_ADD(n,2))" \n"                                             \
                                 : : :BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), CLOTHER_register, ~) "memory"              \
                            );                                                                                      \
                        }                                                                                           \

                     BOOST_PP_REPEAT(MAX_ITERATION_MINUS_ONE, FUNCTION_sub_nbits_nminus1bits, ~)
                     #undef FUNCTION_sub_nbits_nminus1bits

                    } // end namespace detail
             } // end namespace vli
