/*
*Very Large Integer Library, License - Version 1.0 - May 3rd, 2012
*
*Timothee Ewart - University of Geneva, 
*Andreas Hehn - Swiss Federal Institute of technology Zurich.
*
*Permission is hereby granted, free of charge, to any person or organization
*obtaining a copy of the software and accompanying documentation covered by
*this license (the "Software") to use, reproduce, display, distribute,
*execute, and transmit the Software, and to prepare derivative works of the
*Software, and to permit third-parties to whom the Software is furnished to
*do so, all subject to the following:
*
*The copyright notices in the Software and this entire statement, including
*the above license grant, this restriction and the following disclaimer,
*must be included in all copies of the Software, in whole or in part, and
*all derivative works of the Software, unless such copies or derivative
*works are solely in the form of machine-executable object code generated by
*a source language processor.
*
*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
*SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
*FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
*ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
*DEALINGS IN THE SOFTWARE.
*/

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

//g++ -DNUM=1 -E -P -I /opt/boost/include/ main.cpp | sed  "s/n/; \\`echo -e '\n\r      '`/g"
#define MAX_ITERATION 7
#define MAX_ITERATION_MINUS_ONE 6
#define THREE 3
#define FOUR 4
//give the name of the function addition
#define NAME_ADD_NBITS_PLUS_NBITS(n)                 BOOST_PP_CAT(BOOST_PP_CAT(add,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64)),BOOST_PP_CAT(_,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64)))  /* addnx64_nx64 */
//give the name of the multiplication VLI<2*n> = VLI<n>*VLI<n> -  mul2nxx64_nx64_nx64
#define NAME_MUL_TWONBITS_NBITS_NBITS(n) BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_CAT(mul,BOOST_PP_CAT(BOOST_PP_MUL(BOOST_PP_ADD(n,1),2),xx64)),BOOST_PP_CAT(_,BOOST_PP_CAT(BOOST_PP_ADD(n,1),x64))),BOOST_PP_CAT(_,BOOST_PP_CAT(BOOST_PP_ADD(n,1),x64)))  
//give the name of the negation
#define NAME_NEGATE_NBITS(n) BOOST_PP_CAT(neg,BOOST_PP_CAT(BOOST_PP_ADD(n,2),x64))  

//The pp is limited to 256 for arithmetic therefore I calculated intermediate value, close your eyes
//Addition
#define add2x64_2x64 add128_128
#define add3x64_3x64 add192_192
#define add4x64_4x64 add256_256
#define add5x64_5x64 add320_320
#define add6x64_6x64 add384_384
#define add7x64_7x64 add448_448
#define add8x64_8x64 add512_512

//Multiplication
#define mul2xx64_1x64_1x64 mul128_64_64
#define mul4xx64_2x64_2x64 mul256_128_128
#define mul6xx64_3x64_3x64 mul384_192_192
#define mul8xx64_4x64_4x64 mul512_256_256

//negationd
#define neg2x64 neg128
#define neg3x64 neg192
#define neg4x64 neg256
#define neg5x64 neg320
#define neg6x64 neg384
#define neg7x64 neg448
#define neg8x64 neg512

//% not possible with pp
#define pc0 %0
#define pc1 %1
#define pc2 %2
#define pc3 %3
#define pc4 %4
#define pc5 %5
#define pc6 %6
#define pc7 %7
#define pc8 %8
#define pc9 %9
#define pc10 %10
#define pc11 %11
#define pc12 %12
#define pc13 %13
#define pc14 %14
#define pc15 %15
#define pc16 %16
#define pc17 %17
#define pc18 %18
#define pc19 %19

//macro to get the correct name of the register
#define R(n)        BOOST_PP_STRINGIZE(BOOST_PP_CAT(pc,n)) // give register starts from r8 
#define CLOTHER_register_rw(z, n, unused)  BOOST_PP_COMMA_IF(n)"+r"(x[n]) 
#define CLOTHER_register_r(z, n, MAX)  BOOST_PP_COMMA_IF(n)"r"(x[BOOST_PP_ADD(MAX,BOOST_PP_ADD(n,1))]) 

// macro for calculating the indices of the addition
#define I(i,N) BOOST_PP_ADD(i,BOOST_PP_MUL(4,N)) 

// negate for 2CM method, combine with ADC0_register macro
#define NOT_register(z, n, unused)  "not.b32 "R(n)", "R(n)"; \n\t " 
#define ADC0_register(z, n, MAX)    "addc.cc.u32 "R(BOOST_PP_ADD(n,1))", "R(BOOST_PP_ADD(n,1))", "BOOST_PP_STRINGIZE(BOOST_PP_CAT(pc,MAX))"; \n\t " 

//macro for the wrapper

#define FOUR 4 // x,y,z,w
#define VARIABLE0 'x'
#define VARIABLE1 'y'
#define VARIABLE2 'z'
#define VARIABLE3 'w'

//just for fun
//#define QUOTE '
//#define TIM_PP_CHARIZE(ARG) BOOST_PP_CAT(BOOST_PP_CAT(QUOTE,ARG),QUOTE)

#define GET_VAR(z,n,unused) BOOST_PP_COMMA_IF(n) vli::var<BOOST_PP_CAT(VARIABLE,n)> 
#define GET_NULVAR(z,n,unused) BOOST_PP_COMMA_IF(BOOST_PP_SUB(FOUR,n)) vli::no_variable
#define EXPEND_VAR(ARG) BOOST_PP_REPEAT(ARG,GET_VAR,~)  BOOST_PP_REPEAT(BOOST_PP_SUB(FOUR,ARG),GET_NULVAR,~)
//give somehting like  vli::var<'x'> , vli::var<'y'> , vli::no_variable , vli::no_variable
