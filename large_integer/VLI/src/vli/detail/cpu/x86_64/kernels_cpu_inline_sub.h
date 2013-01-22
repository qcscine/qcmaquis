//
//  inline_add.h
//  vli
//
//  Created by Timothée Ewart on 22.01.13.
//
//

#ifndef vli_inline_sub_h
#define vli_inline_sub_h

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

#define EIGHT_MULTIPLY(n)    BOOST_PP_STRINGIZE(BOOST_PP_MUL(n,8))

#define SUBTRACTION( z, n, unused) "movq "EIGHT_MULTIPLY(n)"(%[y]), %[tmp_register] \n\t" \
     BOOST_PP_IF(n,BOOST_PP_STRINGIZE(sbbq),BOOST_PP_STRINGIZE(subq))" %[tmp_register], "EIGHT_MULTIPLY(n)"(%[x]) \n\t" \

#define SUBTRACTION2(z, n, unused) "movq "EIGHT_MULTIPLY(n)"(%[y],%[counter],8), %[tmp_register]\n\t" \
    "sbbq %[tmp_register], "EIGHT_MULTIPLY(n)"(%[x], %[counter], 8)\n\t" \

#define GENERATE_SUBTRACTION(m)  BOOST_PP_REPEAT(BOOST_PP_ADD(m,2), SUBTRACTION, ~)

#define GENERATE_SUBTRACTION2(m)  BOOST_PP_REPEAT(m, SUBTRACTION2, ~)

template<std::size_t n, typename range = void>
struct helper_inline_sub;

/**
 \brief subtraction between two integer<n*128>
 \return void
 \param x boost::uint64_t* pointer of the first entry of the Integer container
 \param y boost::uint64_t* pointer of the first entry of the Integer container
 This operator performs a += between two integer<n*128> The ASM solver is specific and full unroll
 n belongs to the interval 1 to 8, 128-bit <= vli::integer<n*128>  <= 5120-bit
*/
#define FUNCTION_INLINE_sub_nbits_nbits(z, m, unused)                                  \
template<std::size_t n>                                                                \
struct helper_inline_sub<n,typename boost::enable_if_c< n == BOOST_PP_ADD(m,2)>::type>{\
    static void inline_sub(boost::uint64_t* x, boost::uint64_t const* y){              \
        boost::uint64_t tmp_register;                                                  \
            __asm__ __volatile__ (                                                     \
            GENERATE_SUBTRACTION(m)                                                    \
            : [tmp_register] "=&r" (tmp_register)                                      \
            : [x] "r" (x), [y] "r" (y)                                                 \
            : "memory", "cc");                                                         \
    };                                                                                 \
};                                                                                     \

BOOST_PP_REPEAT(8, FUNCTION_INLINE_sub_nbits_nbits, ~) //unroll until 512, maybe test 1024

#undef FUNCTION_INLINE_sub_nbits_nbits


/**
 \brief subtraction between two integer<n*128>, n*128 > 512
 \return void
 \param x boost::uint64_t* pointer of the first entry of the Integer container
 \param y boost::uint64_t* pointer of the first entry of the Integer container
 This operator performs a -= between two integer<n*128> The ASM solver is generic
 and full unroll n belongs to the interval 1 to 8, 512-bit <= vli::integer<n*128>  */
#define FUNCTION_INLINE_sub_nbits_nbits(z, m, unused)                                     \
template<std::size_t n>                                                                   \
struct helper_inline_sub<n,typename boost::enable_if_c<  ((n>8) && (n%8 == m )) >::type>{ \
    static void inline_sub(boost::uint64_t* x, boost::uint64_t const* y){                 \
        boost::uint64_t tmp_register, tmp_register2;                                      \
        __asm__ __volatile__ (                                                           \
            "clc\n\t"                                                                     \
            "1:\n\t"  /* OK I unroll until 8 maybe 16 a day */                            \
            "movq (%[y],%[counter],8), %[tmp_register]\n\t"                               \
            "sbbq %[tmp_register], (%[x], %[counter], 8)\n\t"                             \
                                                                                          \
            "movq 8(%[y],%[counter],8), %[tmp_register]\n\t"                              \
            "sbbq %[tmp_register], 8(%[x], %[counter], 8)\n\t"                            \
                                                                                          \
            "movq 16(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "sbbq %[tmp_register], 16(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 24(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "sbbq %[tmp_register], 24(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 32(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "sbbq %[tmp_register], 32(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 40(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "sbbq %[tmp_register], 40(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 48(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "sbbq %[tmp_register], 48(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "movq 56(%[y],%[counter],8), %[tmp_register]\n\t"                             \
            "sbbq %[tmp_register], 56(%[x], %[counter], 8)\n\t"                           \
                                                                                          \
            "leaq 7(,%[counter],1), %[counter]\n\t" /* lea does not change zero flag so counter += 7 */  \
            "incq %[counter]\n\t" /* change the zero flag */                              \
            "jnz 1b\n\t" /* check if I can leav the loop */                               \
            GENERATE_SUBTRACTION2(m) /* finish the unroll by hand */                      \
        : [tmp_register] "=&r" (tmp_register), "=r" (tmp_register2)                       \
        : [x] "r" (x+(n-m)), [y] "r" (y+(n-m)), [counter] "1" (-(n-m))                    \
        : "memory", "cc");                                                                \
    };                                                                                    \
};                                                                                        \

BOOST_PP_REPEAT(8, FUNCTION_INLINE_sub_nbits_nbits, ~) //unroll until 512, maybe test 1024

#undef FUNCTION_INLINE_sub_nbits_nbits





#undef EIGHT_MULTIPLY
#undef SUBTRACTION
#undef SUBTRACTION2
#undef GENERATE_SUBTRACTION
#undef GENERATE_SUBTRACTION2

#endif
