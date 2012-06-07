#include "vli/utils/macro.h"
// to check :  g++ -DNUM=1 -E -P -I /apps/eiger/boost_1_46_1/include/ -I ../.. vli_number_cpu_function_hooks.hpp | sed  "s/n/;\\`echo -e '\n\r'`/g"  
namespace vli{
    namespace detail{
                     //new functions type : VLI<n*64> *= long int;
                     #define FUNCTION_mul_nbits_64bits(z, n, unused) \
                         void NAME_MUL_NBITS_64BITS(n)(unsigned long int* x, unsigned long int const* y){           \
                         asm(                                                                                       \
                             "movq (%%rsi)          ,%%rax                  \n" /* a0 into rax */                   \
                             "xorq %%rcx            ,%%rcx                  \n" /* rcx to 0 */                      \
                             "cmpq %%rax            ,%%rcx                  \n" /* rax is negative ? */             \
                             "js   "NAME_CONDITIONAL_MUL_NBITS_64BITS(n)"   \n" /* if statements begins */          \
                             "negq %%rax                                    \n" /* negate the number */             \
                             "movq $0x1             ,%%rcx                  \n" /* keep trace for final sign */     \
                             " "NAME_CONDITIONAL_MUL_NBITS_64BITS(n)" :     \n" /* putain de : */                   \
                             "movq %%rax            ,%%rbx                  \n" /* keep a copy of rax/a0 rbx*/      \
                             "mulq "PPS(0,n)"(%%rdi)                        \n" /* lo rax, hi rdx   a0*b0 */        \
                             "movq %%rax            ,%%r8                   \n" /* only one term, write into r8 */  \
                             "movq %%rdx            ,%%r9                   \n" /* hia0b0 into r9 */                \
                             "movq %%rbx            ,%%rax                  \n" /* reload rax */                    \
                              BOOST_PP_REPEAT(n, MUL_register, ~)               /* mul algo */                      \
                             "imulq "PPS(BOOST_PP_ADD(n,1),1)"(%%rdi)       \n" /* a0 * b2, we skip the the hi */   \
                             "addq %%rax            ,"R(BOOST_PP_ADD(n,1))" \n" /* add hi + low */                  \
                             "cmpq $0               ,%%rcx                  \n" /* rcx = 1 we negate */             \
                             "je "NAME_RES_CONDITIONAL_MUL_NBITS_64BITS(n)" \n" /* not equal ZF = 0, negate*/       \
                              BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), NOT_register, ~) /* if for final sign */           \
                             "addq $0x1             ,%%r8                   \n  " /* 2cm add 1 */                   \
                              BOOST_PP_REPEAT(BOOST_PP_ADD(n,1), ADC0_register, ~)/* propagate carry bit */         \
                             " "NAME_RES_CONDITIONAL_MUL_NBITS_64BITS(n)" : \n"   /* end final if */                \
                              BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), SAVE_register, ~)                                  \
                              : : :"rax","rbx","rcx","rdx",BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), CLOTHER_register, ~) "memory"   /* clother register*/      \
                             ); \
                         } \
 
                      BOOST_PP_REPEAT(7, FUNCTION_mul_nbits_64bits, ~) // 7 -> expand until 512 !

                    } // end namespace detail
             } // end namespace vli
