#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/iteration/local.hpp>

#define NAME_MUL_NBITS_64BITS(n)                      BOOST_PP_CAT(BOOST_PP_CAT(mul        ,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64)),BOOST_PP_CAT(_,64)) /* addnx64_64*/
#define NAME_MUL_NBITS_NBITS(n) BOOST_PP_CAT(BOOST_PP_CAT(mul,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64)),BOOST_PP_CAT(_,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64))) /* addnx64_64*/
#define NAME_MUL_TWONBITS_NBITS_NBITS(n)              BOOST_PP_CAT(BOOST_PP_CAT(mul,BOOST_PP_CAT(BOOST_PP_MUL(n,2),x64)),BOOST_PP_CAT(_,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64))) /* addnx64_64*/
#define NAME_CONDITIONAL_MUL_NBITS_64BITS(n)          BOOST_PP_STRINGIZE(BOOST_PP_CAT(BOOST_PP_CAT(_IsNegative,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64)),BOOST_PP_CAT(_,64))) /* _IsNegativenx64_64, for the input sign*/
#define NAME_RES_CONDITIONAL_MUL_NBITS_64BITS(n)      BOOST_PP_CAT(BOOST_PP_CAT(_IsNegativeRes,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64)),BOOST_PP_CAT(_,64)) /* _IsNegativeResnx64_64, for the output sign*/

// boost_pp is limiter to 256 for arithmetic therefore I calculated intermediate value
#define mul2x64_64 mul128_64
#define mul3x64_64 mul192_64
#define mul4x64_64 mul256_64
#define mul5x64_64 mul320_64
#define mul6x64_64 mul384_64
#define mul7x64_64 mul448_64
#define mul8x64_64 mul512_64

#define mul2x64_2x64 mul128_128
#define mul3x64_3x64 mul192_192
#define mul4x64_4x64 mul256_256
#define mul5x64_5x64 mul320_320
#define mul6x64_6x64 mul384_384
#define mul7x64_7x64 mul448_448
#define mul8x64_8x64 mul512_512


#define R(n) BOOST_PP_STRINGIZE(BOOST_PP_CAT(%%r, BOOST_PP_ADD(8,n))) // give register start at r8  

#define PPS(m,n ) BOOST_PP_STRINGIZE( BOOST_PP_MUL(BOOST_PP_MUL(m,n),8)) // m*n*8, 8 because long int
#define PPSr(max,n) BOOST_PP_STRINGIZE( BOOST_PP_MUL(BOOST_PP_ADD(max,n),8)) // m*n*8, 8 because long int

#define LOAD_register(z, n, unused) "movq "PPSr(1,n)"(%%rdi)                 , "R(n)"    \n" /* load 0x??(%%rdi) */     

#define LOADr_register(z, n, unused) "movq "PPSr(1,n)"(%%rsi)                , "Rr(n)"    \n" /* load 0x??(%%rdi), r for reverse */     

#define XOR_register(z, n, unused) "xorq "R(BOOST_PP_ADD(n,2))"           ,"R(BOOST_PP_ADD(n,2))"       \n" /* set up register to 0 */

#define SAVE_register(z, n, unused) "movq "R(n)"                              ,"PPS(AOS,n)"(%%rdi)    \n" /* save 0x??(%%rdi) */     

#define ADC0_register(z, n, unused) "adcq $0x0             ,"R(BOOST_PP_ADD(n,1))"      \n" /* adcq 0 + rdi + CB    */     

#define ADCMUL0_register(z, n, unused) "adcq $0x0             ,"R(BOOST_PP_ADD(n,4))"      \n" /* adcq 0 + rdi + CB    */     

#define  MUL_register(z, n, unused) "mulq "PPS(1,BOOST_PP_ADD(n,1))"(%%rdi)             \n" /* mulq r??*rax */                \
                                    "addq %%rax            ,"R(BOOST_PP_ADD(n,1))"      \n" /* add hia?b? + loa?b? */         \
                                    "movq %%rdx            ,"R(BOOST_PP_ADD(n,2))"      \n" /* save the hi into rcx */        \
                                    "adcq $0               ,"R(BOOST_PP_ADD(n,2))"      \n" /* perhaps carry */               \
                                    "movq %%rbx            ,%%rax                       \n" /* reload rax(a0) from the rbx */ \


#define  MULNN_register(z, n, unused) \
                                      "mulq "PPS(1,n)"(%%rdi)             \n" /* mulq r??*rax */                \
                                      "movq %%rax            ,"R(n)"     \n" /* add hia?b? + loa?b? */         \
                                      "adcq %%rdx            ,"R(BOOST_PP_ADD(n,1))"     \n" /* add hia?b? + loa?b? */         \
                                       BOOST_PP_REPEAT(n, ADCMUL0_register, ~)                                      \

#define NOT_register(z, n, unused)  "notq "R(n)"                                        \n" /* start C2M negate */ 


#define CLOTHER_register(z, n, unused) R(n), /* "r8","r9", ... */

#define Rax(n) BOOST_PP_STRINGIZE(BOOST_PP_CAT(%%r, BOOST_PP_SUB(15,n))) 

 // give register start at r15, r12, .... reverse order  
#define Rdx(n) BOOST_PP_STRINGIZE(BOOST_PP_CAT(%%r, BOOST_PP_SUB(14,n))) // give register start at r15, r12, .... reverse order  
#define Rax(n) BOOST_PP_STRINGIZE(BOOST_PP_CAT(%%r, BOOST_PP_SUB(13,n))) // give register start at r15, r12, .... reverse order  

#define TEST1(z, n, MAX) \
        "movq "Rax(n)", %%rax \n" \
        "movq "Rax(n)", %%rdx \n" \

        BOOST_PP_REPEAT(n, TEST2, unused )                                      \

#define TEST2(z, n, MAX) \

#define BOOST_PP_LOCAL_MACRO(n) \
          -----------------------------------------------------------------" \n"  \
        BOOST_PP_REPEAT(BOOST_PP_ADD(n,2), TEST1, 2) \
   }; \
 

#define BOOST_PP_LOCAL_LIMITS (0, 2)

#include BOOST_PP_LOCAL_ITERATE()



