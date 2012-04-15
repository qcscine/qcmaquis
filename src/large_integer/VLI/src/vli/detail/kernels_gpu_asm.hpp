
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>

// number of iteration for pp add and mul
#define MAX_ITERATION_ADD 2
#define MAX_ITERATION_MUL 6
// macro for calculating the indices of the addition
#define I(i,N) BOOST_PP_ADD(i,BOOST_PP_MUL(6,N)) 

namespace vli{
    namespace detail{
    /* note : ptx is not pure GPU assembly, it is an intermediate,
     * moreover the ptx extented assembly does not support 64 bits integer
     * therefore I have to work on 32 bits, twice more operations for the addition
     * 4 times more for the multuplication
     */    
    inline void add384_384_gpu(unsigned int* x /* shared */, unsigned int const* y /* global */){
    /* this version 60 more ptx lines, so boost pp  
     *  
     * asm( "add.cc.u32   %0 , %0 , %1 ; \n\t" : "+r"(x[0]):"r"(y[0])); 
     * #pragma unroll
     * for(int i=1; i < 11; ++i)
     *     asm( "addc.cc.u32  %0 , %0 , %1 ; \n\t" : "+r"(x[i]):"r"(y[i])); 
     *
     * I have to break up into 2 parts because I can not have more than 30 input/output          
     * load/write operation are done by the compiler (!= ASM x80-86) 
     */
           #define add384_384_384_gpu(w, n, unused) \
               asm( \
                    BOOST_PP_IF(n,"addc.cc.u32 %0, %0, %6; \n\t","add.cc.u32  %0, %0, %6 ; \n\t") /* n=0 no CB,  n!=0, second pass needs CB */ \
                   "addc.cc.u32 %1, %1, %7 ; \n\t" /* x[1] += y[1] + CB                                                                     */ \
                   "addc.cc.u32 %2, %2, %8 ; \n\t" /* x[2] += y[2] + CB                                                                     */ \
                   "addc.cc.u32 %3, %3, %9 ; \n\t" /* x[3] += y[3] + CB                                                                     */ \
                   "addc.cc.u32 %4, %4, %10; \n\t" /* x[4] += y[4] + CB                                                                     */ \
                   "addc.cc.u32 %5, %5, %11; \n\t" /* x[5] += y[5] + CB                                                                     */ \
                   :"+r"(x[I(0,n)]),"+r"(x[I(1,n)]),"+r"(x[I(2,n)]),                                                                           \
                    "+r"(x[I(3,n)]),"+r"(x[I(4,n)]),"+r"(x[I(5,n)])                                                                            \
                   :"r"(y[I(0,n)]),"r"(y[I(1,n)]),"r"(y[I(2,n)]),                                                                              \
                    "r"(y[I(3,n)]),"r"(y[I(4,n)]),"r"(y[I(5,n)])                                                                               \
                  );
           BOOST_PP_REPEAT(MAX_ITERATION_ADD, add384_384_384_gpu, ~)

           #undef add384_384_384_gpu
    } //end add384_384_gpu 

    inline void mul384_384_gpu(unsigned int* x/* res local*/, unsigned int const* y /* shared */, unsigned int const* z /* shared */){
    /*                   y[5] y[4] y[3] y[2] y[1] y[0]  %13 %12 %11 %10 %9 %8  
     *                 X                          z[i], %7  i = 0 in my example
     * _________________________________________________
     * x[i+6] x[i+5] x[i+4] x[i+3] x[i+2] x[i+1] x[i+0], %6 %5 %4 %3 %2 %1 %0
     * I multiply first low part and then high part, I avoid all carry bit propagation pbs
     * The pp_if CB only possible when n!=0 because z is init to 0 by default
     */
           #define mul384_192_192_gpu(w, n, unused) \
               asm( \
                  "mad.lo.cc.u32  %0, %8,  %7, %0; \n\t" /* c[i]   = a[0] * b[i] (low)  + c[i] (c[i]=0 for i=0) may generate carry bit (CB) */ \
                  "madc.lo.cc.u32 %1, %9,  %7, %1; \n\t" /* c[i+1] = a[1] * b[i] (low)  + c[i+1] + CB                                       */ \
                  "madc.lo.cc.u32 %2, %10, %7, %2; \n\t" /* c[i+2] = a[2] * b[i] (low)  + c[i+1] + CB                                       */ \
                  "madc.lo.cc.u32 %3, %11, %7, %3; \n\t" /* c[i+3] = a[3] * b[i] (low)  + c[i+3] + CB                                       */ \
                  "madc.lo.cc.u32 %4, %12, %7, %4; \n\t" /* c[i+4] = a[4] * b[i] (low)  + c[i+4] + CB                                       */ \
                  "madc.lo.cc.u32 %5, %13, %7, %5; \n\t" /* c[i+5] = a[5] * b[i] (low)  + c[i+5] + CB                                       */ \
                  BOOST_PP_IF(n,"addc.cc.u32 %6, 0, 0; \n\t",/* no extention n=0 */) /* c[i+6] += CB, n = 0 CB impossible                   */ \
                  "mad.hi.cc.u32  %1, %8,  %7, %1; \n\t" /* c[i+1] = a[0] * b[i] (high) + c[i+1] + CB (c[i]=0 for i=0)                      */ \
                  "madc.hi.cc.u32 %2, %9,  %7, %2; \n\t" /* c[i+2] = a[1] * b[i] (high) + c[i+2] + CB                                       */ \
                  "madc.hi.cc.u32 %3, %10, %7, %3; \n\t" /* c[i+3] = a[2] * b[i] (high) + c[i+3] + CB                                       */ \
                  "madc.hi.cc.u32 %4, %11, %7, %4; \n\t" /* c[i+4] = a[3] * b[i] (high) + c[i+4] + CB                                       */ \
                  "madc.hi.cc.u32 %5, %12, %7, %5; \n\t" /* c[i+5] = a[4] * b[i] (high) + c[i+5] + CB                                       */ \
                  "madc.hi.cc.u32 %6, %13, %7, %6; \n\t" /* c[i+6] = a[5] * b[i] (high) + c[i+6]                                            */ \
                  :"+r"(x[BOOST_PP_ADD(0,n)]),"+r"(x[BOOST_PP_ADD(1,n)]),"+r"(x[BOOST_PP_ADD(2,n)]),                                           \
                   "+r"(x[BOOST_PP_ADD(3,n)]),"+r"(x[BOOST_PP_ADD(4,n)]),"+r"(x[BOOST_PP_ADD(5,n)]),"+r"(x[BOOST_PP_ADD(6,n)])                 \
                  :"r"(z[n]),"r"(y[0]),"r"(y[1]),"r"(y[2]),"r"(y[3]),"r"(y[4]),"r"(y[5])                                                       \
                 );                                                                                                                            \

           BOOST_PP_REPEAT(MAX_ITERATION_MUL, mul384_192_192_gpu, ~)
           #undef mul384_192_192_gpu
   } // end mul384_384_gpu

   inline void muladd384_384_gpu(unsigned int* x/* res local*/, unsigned int const* y /* shared */, unsigned int const* z /* shared */){

   }

   } //end namespace detail
} //end namespace vli
