//
//  inline_add.h
//  vli
//
//  Created by Timothée Ewart on 22.01.13.
//
//

#ifndef vli_inline_add_h
#define vli_inline_add_h

//SFINAE  http://www.boost.org/doc/libs/1_52_0/libs/utility/enable_if.html


//Address operand syntax, I always forget
//
//There are up to 4 parameters of an address operand that are presented in the syntax displacement(base register, offset register, scalar multiplier). This is equivalent to [base register + displacement + offset register * scalar multiplier] in Intel syntax. Either or both of the numeric, and either of the register parameters may be omitted:
//movl    -4(%ebp, %edx, 4), %eax  # Full example: load *(ebp - 4 + (edx * 4)) into eax
//movl    -4(%ebp), %eax           # Typical example: load a stack variable into eax
//movl    (%ecx), %edx             # No offset: copy the target of a pointer into a register
//leal    8(,%eax,4), %eax         # Arithmetic: multiply eax by 4 and add 8
//leal    (%eax,%eax,2), %eax      # Arithmetic: multiply eax by 2 and add eax (i.e. multiply by 3)

// "movq (%[y],%[counter],8), %[tmp_register]\n\t" load *(y - 0 + counter*8)


#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/stringize.hpp>

#define VLI_EIGHT_MULTIPLY(n)    BOOST_PP_STRINGIZE(BOOST_PP_MUL(n,8))

#define VLI_ADDITION( z, n, unused) "movq "VLI_EIGHT_MULTIPLY(n)"(%[y]), %[tmp_register] \n\t" \
     BOOST_PP_IF(n,BOOST_PP_STRINGIZE(adcq),BOOST_PP_STRINGIZE(addq))" %[tmp_register], "VLI_EIGHT_MULTIPLY(n)"(%[x]) \n\t" \

#define VLI_ADDITION2(z, n, unused) "movq "VLI_EIGHT_MULTIPLY(n)"(%[y],%[counter],8), %[tmp_register]\n\t" \
                                "adcq %[tmp_register], "VLI_EIGHT_MULTIPLY(n)"(%[x], %[counter], 8)\n\t" \

#define VLI_ADDITION3(z, n, unused) \
     BOOST_PP_IF(n,BOOST_PP_STRINGIZE(adcq %[constante2]),BOOST_PP_STRINGIZE(addq %[tmp_register]))", "VLI_EIGHT_MULTIPLY(n)"(%[x]) \n\t"

#define VLI_ADDITION4(z, n, unused) "adcq %[constante2], "EIGHT_MULTIPLY(n)"(%[x], %[counter], 8)\n\t" \

#define VLI_GENERATE_ADDITION(m)  BOOST_PP_REPEAT(BOOST_PP_ADD(m,2), VLI_ADDITION, ~)

#define VLI_GENERATE_ADDITION2(m)  BOOST_PP_REPEAT(m, VLI_ADDITION2, ~)

#define VLI_GENERATE_ADDITION3(m)  BOOST_PP_REPEAT(BOOST_PP_ADD(m,2), VLI_ADDITION3, ~)

#define VLI_GENERATE_ADDITION4(m)  BOOST_PP_REPEAT(m, VLI_ADDITION4, ~)

template<std::size_t NumWords, typename range = void>
struct helper_inline_add;

template<std::size_t NumWords, typename range = void>
struct helper_inline_add_constante;

/**
 \brief addition between two integer<n*128>, n*128 <= 512
 \return void
 \param x boost::uint64_t* pointer of the first entry of the Integer container
 \param y boost::uint64_t* pointer of the first entry of the Integer container
 This operator performs a += between two integer<n*128> The ASM solver is specific 
 and full unroll n belongs to the interval 1 to 8, 128-bit <= vli::integer<n*128>  <= 512-bit
*/
#define FUNCTION_INLINE_add_nbits_nbits(z, m, unused)                                  \
template<std::size_t NumWords>                                                         \
struct helper_inline_add<NumWords,typename boost::enable_if_c< NumWords == BOOST_PP_ADD(m,2)>::type>{\
    static void inline_add(boost::uint64_t* x, boost::uint64_t const* y){              \
        boost::uint64_t tmp_register;                                                  \
        __asm__ __volatile__ (                                                         \
            VLI_GENERATE_ADDITION(m)                                                   \
        : [tmp_register] "=&r" (tmp_register)                                          \
        : [x] "r" (x), [y] "r" (y)                                                     \
        : "memory", "cc");                                                             \
    };                                                                                 \
                                                                                       \
    static void inline_add(boost::uint64_t* x, boost::int64_t const y){                \
        boost::uint64_t tmp_register;                                                  \
        boost::uint64_t tmp_register2(y >> 63);                                        \
        __asm__ __volatile__ (                                                         \
            "movq  %[constante],   %[tmp_register] \n\t"                               \
            VLI_GENERATE_ADDITION3(m)                                                  \
        : [tmp_register] "=&r" (tmp_register)                                          \
        : [x] "r" (x), [constante] "r" (y), [constante2] "r" (tmp_register2)           \
        : "memory", "cc");                                                             \
    };                                                                                 \
};                                                                                     \

BOOST_PP_REPEAT(8, FUNCTION_INLINE_add_nbits_nbits, ~) //unroll until 512, maybe test 1024

#undef FUNCTION_INLINE_add_nbits_nbits


/**
 \brief addition between two integer<n*128>, n*128 > 512
 \return void
 \param x boost::uint64_t* pointer of the first entry of the Integer container
 \param y boost::uint64_t* pointer of the first entry of the Integer container
 This operator performs a += between two integer<n*128> The ASM solver is specific
 and full unroll n belongs to the interval 1 to 8, 512-bit <= vli::integer<n*128>  */
#define FUNCTION_INLINE_add_nbits_nbits(z, m, unused)                                     \
template<std::size_t NumWords>                                                                   \
struct helper_inline_add<NumWords,typename boost::enable_if_c<  ((NumWords>8) && (NumWords%8 == m )) >::type>{ \
    static void inline_add(boost::uint64_t* x, boost::uint64_t const* y){                 \
        boost::uint64_t tmp_register, tmp_register2;                                      \
        __asm__ __volatile__ (                                                            \
            "clc\n\t"                                                                     \
            "1:\n\t"  /* OK I unroll until 8 maybe 16 a day */                            \
            "movq (%[y],%[counter],8), %[tmp_register]\n\t"                               \
            "adcq %[tmp_register], (%[x], %[counter], 8)\n\t"                             \
                                                                                          \
            "movq 8(%[y],%[counter],8), %[tmp_register]\n\t"                              \
            "adcq %[tmp_register], 8(%[x], %[counter], 8)\n\t"                            \
                                                                                          \
            "movq 16(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "adcq %[tmp_register], 16(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 24(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "adcq %[tmp_register], 24(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 32(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "adcq %[tmp_register], 32(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 40(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "adcq %[tmp_register], 40(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 48(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "adcq %[tmp_register], 48(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 56(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "adcq %[tmp_register], 56(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "leaq 7(,%[counter],1), %[counter]\n\t" /* lea does not change zero flag so counter += 7 */  \
            "incq %[counter]\n\t" /* change the zero flag */                              \
            "jnz 1b\n\t" /* check if I can leav the loop */                               \
            VLI_GENERATE_ADDITION2(m) /* finish the unroll by hand */                     \
        : [tmp_register] "=&r" (tmp_register), "=r" (tmp_register2)                       \
        : [x] "r" (x+(NumWords-m)), [y] "r" (y+(NumWords-m)), [counter] "1" (-(NumWords-m))\
        : "memory", "cc");                                                                \
    };                                                                                    \
};                                                                                        \

BOOST_PP_REPEAT(8, FUNCTION_INLINE_add_nbits_nbits, ~) //unroll until 512, maybe test 1024

#undef FUNCTION_INLINE_add_nbits_nbits

/*

template<std::size_t NumWords>
struct helper_inline_add_constante<NumWords,typename boost::enable_if_c<2>::type>{
    static void inline_add_constante(boost::uint64_t* x, boost::int64_t const y){
        boost::uint64_t tmp_register;
        __asm__ __volatile__ (
                 "movq  %[constante] ,   %[tmp_register] \n\t"
                 "sarq   $63           , %[constante]    \n\t"
                 "addq  %[tmp_register], (%[x]) \n\t"
                 "adcq  %[constante],  8(%[x]) \n\t"
        : [tmp_register] "=&r" (tmp_register)
        : [x] "r" (x), [constante] "r" (y)
        : "memory", "cc");
    };
};

*/



/**
 \brief addition between two integer<n*512>
 \return void
 \param x boost::uint64_t* pointer of the first entry of the Integer container
 \param y boost::uint64_t* pointer of the first entry of the Integer container
 This operator performs a += between two integer<n*512> The ASM solver is generic n != 0
 */

/*
template<std::size_t n>
struct helper_inline_add<n,typename boost::enable_if_c< ((n>8) && (n%8 == 0 )) >::type>{
    static void inline_add(boost::uint64_t* x, boost::uint64_t const* y){
       boost::uint64_t tmp_register, tmp_register2;
        __asm__ __volatile__ (
                              "clc\n\t"
                              "1:\n\t"
                              "movq (%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], (%[x], %[counter], 8)\n\t"

                              "movq 8(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 8(%[x], %[counter], 8)\n\t"
                              
                              "movq 16(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 16(%[x], %[counter], 8)\n\t"

                              "movq 24(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 24(%[x], %[counter], 8)\n\t"
                              
                              "movq 32(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 32(%[x], %[counter], 8)\n\t"
                              
                              "movq 40(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 40(%[x], %[counter], 8)\n\t"
                              
                              "movq 48(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 48(%[x], %[counter], 8)\n\t"
                              
                              "movq 56(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 56(%[x], %[counter], 8)\n\t"
                              
                              "leaq 7(,%[counter],1), %[counter]\n\t" // counter*1+7, this does not change zero flag so I can not make +8, I have to make +7 and inc
                              "incq %[counter]\n\t" // change the zero flag
                              "jnz 1b\n\t" // check if I can leav the loop
                              : [tmp_register] "=&r" (tmp_register), "=r" (tmp_register2)
                              : [x] "r" (x+n), [y] "r" (y+n), [counter] "1" (-n) // ^_^ incq not decq so ...
                              : "memory", "cc");
    };
};

template<std::size_t n>
struct helper_inline_add<n,typename boost::enable_if_c< ((n>8) && (n%8 == 1 )) >::type>{
    static void inline_add(boost::uint64_t* x, boost::uint64_t const* y){
        boost::uint64_t tmp_register, tmp_register2;

        __asm__ __volatile__ (
                              "clc\n\t"
                              "1:\n\t"
                              "movq (%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], (%[x], %[counter], 8)\n\t"
                              
                              "movq 8(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 8(%[x], %[counter], 8)\n\t"
                              
                              "movq 16(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 16(%[x], %[counter], 8)\n\t"
                              
                              "movq 24(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 24(%[x], %[counter], 8)\n\t"
                              
                              "movq 32(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 32(%[x], %[counter], 8)\n\t"
                              
                              "movq 40(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 40(%[x], %[counter], 8)\n\t"
                              
                              "movq 48(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 48(%[x], %[counter], 8)\n\t"
                              
                              "movq 56(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 56(%[x], %[counter], 8)\n\t"
                              
                              "leaq 7(,%[counter],1), %[counter]\n\t"
                              "incq %[counter]\n\t"
                              "jnz 1b\n\t"

                              "movq (%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], (%[x], %[counter], 8)\n\t"
                              
                              : [tmp_register] "=&r" (tmp_register), "=r" (tmp_register2)
                              : [x] "r" (x+(n-1)), [y] "r" (y+(n-1)), [counter] "1" (-(n-1))
                              : "memory", "cc");
    };
};

template<std::size_t n>
struct helper_inline_add<n,typename boost::enable_if_c< ((n>8) && (n%8 == 2 )) >::type>{
    static void inline_add(boost::uint64_t* x, boost::uint64_t const* y){
        boost::uint64_t tmp_register, tmp_register2;
        
        __asm__ __volatile__ (
                              "clc\n\t"
                              "1:\n\t"
                              "movq (%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], (%[x], %[counter], 8)\n\t"
                              
                              "movq 8(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 8(%[x], %[counter], 8)\n\t"
                              
                              "movq 16(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 16(%[x], %[counter], 8)\n\t"
                              
                              "movq 24(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 24(%[x], %[counter], 8)\n\t"
                              
                              "movq 32(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 32(%[x], %[counter], 8)\n\t"
                              
                              "movq 40(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 40(%[x], %[counter], 8)\n\t"
                              
                              "movq 48(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 48(%[x], %[counter], 8)\n\t"
                              
                              "movq 56(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 56(%[x], %[counter], 8)\n\t"
                              
                              "leaq 7(,%[counter],1), %[counter]\n\t"
                              "incq %[counter]\n\t"
                              "jnz 1b\n\t"

                              "movq (%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], (%[x], %[counter], 8)\n\t"
                              
                              "movq 8(%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], 8(%[x], %[counter], 8)\n\t"
                              
                              : [tmp_register] "=&r" (tmp_register), "=r" (tmp_register2)
                              : [x] "r" (x+(n-2)), [y] "r" (y+(n-2)), [counter] "1" (-(n-2)) 
                              : "memory", "cc");
    };
};

*/

/*
template<std::size_t n>
struct helper_inline_add<n,typename boost::enable_if_c< ((n>4) && (n%4 == 1)) >::type>{
    static void inline_add(boost::uint64_t* x, boost::uint64_t const* y){
        std::cout << " je suis la " <<((n/4)*4) << std::endl;
        boost::uint64_t tmp_register, tmp_register2;
        __asm__ __volatile__ (
                              "clc\n\t"
                              "1:\n\t"
                              "movq (%[y],%[counter],8), %[tmp_register]\n\t"
                              "adcq %[tmp_register], (%[x], %[counter], 8)\n\t"
                              "incq %[counter]\n\t"

                              "jnz 1b\n\t"
                              
                              
    
                              
                              : [tmp_register] "=&r" (tmp_register), "=r" (tmp_register2)
                              : [x] "r" (x+n), [y] "r" (y+n), [counter] "1" (-n) // ^_^ incq not decq so ...
                              : "memory", "cc");
    };
};

template<std::size_t n>
struct helper_inline_add<n,typename boost::enable_if_c< ((n>4) && (n%4 == 2 )) >::type>{
    static void inline_add(boost::uint64_t* x, boost::uint64_t const* y){
        std::cout << " JE SUIS LA q*256 + 2" << std::endl;
    };
};

template<std::size_t n>
struct helper_inline_add<n,typename boost::enable_if_c< ((n>4) && (n%4 == 3 )) >::type>{
    static void inline_add(boost::uint64_t* x, boost::uint64_t const* y){
        std::cout << " JE SUIS LA q*256 + 3" << std::endl;
    };
};

*/


/*
 template<long n>
 void inline_add(boost::uint64_t* x, boost::uint64_t const* y) {
 boost::uint64_t tmp_register, tmp_register2;
 __asm__ __volatile__ (
 "clc\n\t"
 "1:\n\t"
 "movq (%[y],%[counter],8), %[tmp_register]\n\t"
 "adcq %[tmp_register], (%[x], %[counter], 8)\n\t"
 "incq %[counter]\n\t"
 "jnz 1b\n\t"
 : [tmp_register] "=&r" (tmp_register), "=r" (tmp_register2)
 : [x] "r" (x), [y] "r" (y), [counter] "1" (n)
 : "memory", "cc");
 }
 */

#undef EIGHT_MULTIPLY
#undef ADDITION
#undef ADDITION2
#undef ADDITION3
#undef ADDITION4

#undef GENERATE_ADDITION
#undef GENERATE_ADDITION2
#undef GENERATE_ADDITION3
#undef GENERATE_ADDITION4

#endif
