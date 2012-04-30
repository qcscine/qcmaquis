//
//  kernels_cpu.cpp
//  VLI_ASM
//
//  Created by Tim Ewart on 22.03.12.
//  Copyright (c) 2012 __Université de Genève__. All rights reserved.
//
#include "kernels_cpu_asm.h"
#include "vli/utils/macro.h"
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/stringize.hpp>

namespace vli{
    namespace detail{
        
// c- this assembly code support both data layout SoA and AoS
// clear the syntax
//#define PPS(m,n) BOOST_PP_STRINGIZE( BOOST_PP_MUL(BOOST_PP_MUL(m,n),8)) // m*n*8, 8 because long int
// PPS(1,n) = 0x08 hex = dec 8
// PPS(2,n) = 0x10 hex = dec 16
// PPS(3,n) = 0x18 hex = dec 24
// PPS(4,n) = 0x20 hex = dec 32
// PPS(5,n) = 0x28 hex = dec 40

// TO DO + - and basic * could boost pp to until 512, to generate all variations

// Note : for the sign
// 192 * 64  = 192 : necessitate to adapt (if), use for VLI *= long , VLI * long
// 384 * 64  = 384 : necessitate to adapt (if), use for print
// 192 * 192 = 192 : natively compatible VLI *= VLI, VLI*VLI
// 192 * 192 = 384 : necessitate to adapt : mul and muladd
        
        
// Note, the ASM X86-64 allows a different managment of the stack
// I am presently using the red zone there I do not have to allocate 
// the stack pointer (the frame pointer is removed under x86-64)
// the red zone is a zone of 128 bytes, enough for my arithmetic        
//


void mul128_64_64(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */, unsigned long int const* z/* %%rdx -> rbx */){
   asm( 
       "movq (%%rsi)          ,%%rax             \n" /* a0 into rax */                   
       "movq %%rdx            ,%%rcx             \n" /* save a0-rcx faster than stack */ 
       "imulq (%%rcx)                            \n" /* lo rax, hi rdx   a0*b0 */        
       "movq %%rax            ,(%%rdi)           \n" /* lo part save */
       "movq %%rdx            ,8(%%rdi)          \n" /* hi part sve */
       : : :"rax","rdx","rcx","memory"
       );
}

void mul256_128_128(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */, unsigned long int const* z/* %%rdx -> rbx */){
   asm( 
/*-01*/ "subq $0x20            ,%%rsp             \n" /* create stack frame */            
/*00*/  "movq %%rdx            ,%%rbx             \n" /* rdx uses by mul             */   
/*01*/  "xorq %%r10            ,%%r10             \n" /* r10 = 0 due to carry effect */   
/*02*/  "xorq %%r11            ,%%r11             \n" /* r11 = 0 due to carry effect */   
/*03*/  "xorq %%r14            ,%%r14             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/*03*/  "xorq %%r15            ,%%r15             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/* negate a if negative and store into stack, in reverse order due to universal access */ 
/*04*/  "movq 8(%%rsi)         ,%%rax             \n" /* load a3 into r8, for the sign */ 
/*05*/  "cmpq $0               ,%%rax             \n" /* test a is negative(sign on a3)*/ 
/*06*/  "jns _Negativea_256_128_                  \n" /* number is negative, we negate */ 
/*07*/  "movq (%%rsi)          ,%%r8              \n" /* load a0 */                       
/*09*/  "notq %%r8                                \n" /* C2M, ~a0 */                      
/*11*/  "notq %%rax                               \n" /* C2M, ~a1 */                      
/*12*/  "addq $0x1             ,%%r8              \n" /* C2M, ~a0+1 */                    
/*14*/  "adcq $0x0             ,%%rax             \n" /* C2M, ~a1+CB */                   
/*16*/  "movq %%r8             ,-0x10(%%rsp)      \n" /* a0 into the stack -16 rsp */     
/*17*/  "movq %%rax            ,-0x08(%%rsp)      \n" /* a1 into the stack -8 rsp */      
/*18*/  "leaq  -0x10(%%rsp)    ,%%rsi             \n" /* rsi points to stack a0 > 0 */    
/*19*/  "movq $1               ,%%r14             \n" /* r14 = 0 it is the sign 0+ 1-*/   
/*20*/  "_Negativea_256_128_ :                    \n" /* end if structure */              
/* negate a if negative and store into stack, in reverse order due to universal access */ 
/*04*/  "movq 8(%%rbx)         ,%%rax             \n" /* load a3 into r8, for the sign */ 
/*05*/  "cmpq $0               ,%%rax             \n" /* test a is negative(sign on a3)*/ 
/*06*/  "jns _Negativeb_256_128_                  \n" /* number is negative, we negate */ 
/*07*/  "movq (%%rbx)          ,%%r8              \n" /* load b0 */                       
/*09*/  "notq %%r8                                \n" /* C2M, ~b0 */                      
/*11*/  "notq %%rax                               \n" /* C2M, ~b1 */                      
/*12*/  "addq $0x1             ,%%r8              \n" /* C2M, ~b0+1 */                    
/*14*/  "adcq $0x0             ,%%rax             \n" /* C2M, ~b1+CB */                   
/*16*/  "movq %%r8             ,-0x20(%%rsp)      \n" /* b0 into the stack -16 rsp */     
/*17*/  "movq %%rax            ,-0x18(%%rsp)      \n" /* b1 into the stack -8 rsp */      
/*18*/  "leaq  -0x20(%%rsp)    ,%%rsi             \n" /* rsi points to stack b0 > 0 */    
/*19*/  "movq $1               ,%%r15             \n" /* r15 = 0 it is the sign 0+ 1-*/   
/*20*/  "_Negativeb_256_128_ :                    \n" /* end if structure */              
/* --------------------------- a0 * b0, a0 * b1 start ------------------------*/ 
/*39*/  "movq (%%rsi)          ,%%rax             \n" /* a0 into rax */                   
/*40*/  "movq %%rax            ,%%rcx             \n" /* save a0-rcx faster than stack */ 
/*41*/  "mulq (%%rbx)                             \n" /* lo rax, hi rdx   a0*b0 */        
/*42*/  "movq %%rax            ,%%r8              \n" /* only one term, write into c0 */  
/*43*/  "movq %%rdx            ,%%r9              \n" /* a0b0hi into r8 */                
/*44*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a0) from the stack */ 
/*45*/  "mulq 8(%%rbx)                            \n" /* a0 * b1 */                       
/*46*/  "addq %%rax            ,%%r9              \n" /* add a0b0hi + a0b1lo */           
/*47*/  "adcq %%rdx            ,%%r10             \n" /* save the a0b1hi into r9 */       
/* --------------------------- a0 * b0, a0 * b1 end --------------------------*/ 
/* --------------------------- a1 * b0, a1 * b1 start ------------------------*/ 
/*52*/  "movq 8(%%rsi)         ,%%rax             \n" /* a1 into rax */                   
/*53*/  "movq %%rax            ,%%rcx             \n" /* save a0-rcx faster than stack */ 
/*54*/  "mulq (%%rbx)                             \n" /* a1 * b0 */                       
/*55*/  "addq %%rax            ,%%r9              \n" /* l46 + a1b0lo */                  
/*56*/  "adcq %%rdx            ,%%r10             \n" /* l47 + a1b0hi + c */              
/*57*/  "adcq $0               ,%%r11             \n" /* possible carry, for 192 one adcq */ 
/*58*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a1) from the stack */ 
/*63*/  "mulq 8(%%rbx)                            \n" /* a1*b1 */                         
/*64*/  "addq %%rax            ,%%r10             \n" /* a1b2lo to r9 */                  
/*65*/  "adcq %%rdx            ,%%r11             \n" /* a1b2hi + c  */                   
/* --------------------------- a1 * b0, a1 * b1 end --------------------------*/ 
/*81*/  "xorq %%r14            ,%%r15             \n"                                     
/*82*/  "cmpq $0               ,%%r15             \n" /* r15 = 1 we negate */             
/*83*/  "je _IsNegativeResult_256_128_            \n" /* not equal ZF = 0, negate*/       
/*84*/  "notq %%r8                                \n" /* start2ComplementMethod negate */ 
/*85*/  "notq %%r9                                \n" /* 2CM negate */                    
/*86*/  "notq %%r10                               \n" /* 2CM negate */                    
/*87*/  "notq %%r11                               \n" /* 2CM negate */                    
/*90*/  "addq $0x1             ,%%r8              \n" /* 2CM add 1 */                     
/*91*/  "adcq $0x0             ,%%r9              \n" /* 2CM propagate CB */              
/*92*/  "adcq $0x0             ,%%r10             \n" /* 2CM propagate CB */              
/*93*/  "adcq $0x0             ,%%r11             \n" /* 2CM propagate CB */              
/*96*/  "_IsNegativeResult_256_128_ :             \n" /* end if negative result */        
/*97*/  "movq %%r8             ,(%%rdi)           \n" /* r8 -> c1 */                      
/*98*/  "movq %%r9             ,8(%%rdi)          \n" /* r9 -> c1 */                      
/*99*/  "movq %%r10            ,16(%%rdi)         \n" /* r10 -> c2 */                     
/*100*/ "movq %%r11            ,24(%%rdi)         \n" /* r11 -> c3 */                     
/*103*/ "addq $0x20            ,%%rsp             \n" /* destroy stack frame */           
        : : : "rax","rbx","rcx","rdx","r8","r9","r10","r11","r14","r15","memory"   
       );
}


//void mul512_256_256(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */, unsigned long int const* z/* %%rdx -> rbx */){
//   asm( 
//
//        : : : "rax","rbx","rcx","rdx","r8","r9","r10","r11","r14","r15","memory"   
//   ); 
//}

#define HELPER_ASM_MUL384_192_192(n) \
void mul384_192_192(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */, unsigned long int const* z/* %%rdx -> rbx */){ \
   asm( \
/*-01*/ "subq $0x30            ,%%rsp             \n" /* create stack frame */            \
/*00*/  "movq %%rdx            ,%%rbx             \n" /* rdx uses by mul             */   \
/*01*/  "xorq %%r11            ,%%r11             \n" /* r10 = 0 due to carry effect */   \
/*02*/  "xorq %%r12            ,%%r12             \n" /* r11 = 0 due to carry effect */   \
/*03*/  "xorq %%r13            ,%%r13             \n" /* r12 = 0 due to carry effect */   \
/*03*/  "xorq %%r14            ,%%r14             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/*03*/  "xorq %%r15            ,%%r15             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/* negate a if negative and store into stack, in reverse order due to universal access */ \
/*04*/  "movq "PPS(2,n)"(%%rsi),%%rax             \n" /* load a3 into r8, for the sign */ \
/*05*/  "cmpq $0               ,%%rax             \n" /* test a is negative(sign on a3)*/ \
/*06*/  "jns _Negativea_384_192_                  \n" /* number is negative, we negate */ \
/*07*/  "movq (%%rsi)          ,%%r8              \n" /* load a0 */                       \
/*08*/  "movq "PPS(1,n)"(%%rsi),%%r9              \n" /* load a1 */                       \
/*09*/  "notq %%r8                                \n" /* C2M, ~a0 */                      \
/*10*/  "notq %%r9                                \n" /* C2M, ~a1 */                      \
/*11*/  "notq %%rax                               \n" /* C2M, ~a2 */                      \
/*12*/  "addq $0x1             ,%%r8              \n" /* C2M, ~a0+1 */                    \
/*13*/  "adcq $0x0             ,%%r9              \n" /* C2M, ~a1+CB */                   \
/*14*/  "adcq $0x0             ,%%rax             \n" /* C2M, ~a2+CB */                   \
/*15*/  "movq %%r8             ,-0x18(%%rsp)      \n" /* a0 into the stack -24 rsp */     \
/*16*/  "movq %%r9             ,-0x10(%%rsp)      \n" /* a1 into the stack -16 rsp */     \
/*17*/  "movq %%rax            ,-0x08(%%rsp)      \n" /* a2 into the stack -8 rsp */      \
/*18*/  "leaq  -0x18(%%rsp)    ,%%rsi             \n" /* rsi points to stack a0 > 0 */    \
/*19*/  "movq $1               ,%%r14             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/*20*/  "_Negativea_384_192_ :                    \n" /* end if structure */              \
/* negate a if negative and store into stack, in reverse order due to universal access */ \
/*21*/  "movq "PPS(2,n)"(%%rbx),%%rax             \n" /* load a3 into r8, for the sign */ \
/*22*/  "cmpq $0               ,%%rax             \n" /* test a is negative(sign on a3)*/ \
/*23*/  "jns _Negativeb_384_192_                  \n" /* number is negative, we negate */ \
/*24*/  "movq (%%rbx)          ,%%r8              \n" /* load b0 */                       \
/*25*/  "movq "PPS(1,n)"(%%rbx),%%r9              \n" /* load b1 */                       \
/*26*/  "notq %%r8                                \n" /* C2M, ~b0 */                      \
/*27*/  "notq %%r9                                \n" /* C2M, ~b1 */                      \
/*28*/  "notq %%rax                               \n" /* C2M, ~b2 */                      \
/*29*/  "addq $0x1             ,%%r8              \n" /* C2M, ~b0+1 */                    \
/*30*/  "adcq $0x0             ,%%r9              \n" /* C2M, ~b1+CB */                   \
/*31*/  "adcq $0x0             ,%%rax             \n" /* C2M, ~b2+CB */                   \
/*32*/  "movq %%r8             ,-0x30(%%rsp)      \n" /* b0 into the stack -48 rsp */     \
/*33*/  "movq %%r9             ,-0x28(%%rsp)      \n" /* b1 into the stack -40 rsp */     \
/*34*/  "movq %%rax            ,-0x20(%%rsp)      \n" /* b2 into the stack -32 rsp */     \
/*35*/  "leaq  -0x30(%%rsp)    ,%%rbx             \n" /* rsi points to stack b0 > 0 */    \
/*36*/  "movq $1               ,%%r15             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/*37*/  "_Negativeb_384_192_ :                    \n" /* end if structure */              \
/*38*/  "xorq %%r10            ,%%r10             \n" /* r9 = 0  due to carry effect */   \
/* --------------------------- a0 * b0, a0 * b1, a0 * b2 start ------------------------*/ \
/*39*/  "movq (%%rsi)          ,%%rax             \n" /* a0 into rax */                   \
/*40*/  "movq %%rax            ,%%rcx             \n" /* save a0-rcx faster than stack */ \
/*41*/  "mulq (%%rbx)                             \n" /* lo rax, hi rdx   a0*b0 */        \
/*42*/  "movq %%rax            ,%%r8              \n" /* only one term, write into c0 */  \
/*43*/  "movq %%rdx            ,%%r9              \n" /* a0b0hi into r8 */                \
/*44*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a0) from the stack */ \
/*45*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a0 * b1 */                       \
/*46*/  "addq %%rax            ,%%r9              \n" /* add a0b0hi + a0b1lo */           \
/*47*/  "adcq %%rdx            ,%%r10             \n" /* save the a0b1hi into r9 */       \
/*48*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a0) from the stack */ \
/*49*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a0 * b2 */                       \
/*50*/  "addq %%rax            ,%%r10             \n" /* add l11 + a0b2lo + c */          \
/*51*/  "adcq %%rdx            ,%%r11             \n" /* add l01 + a0b2hi + c, end */     \
/* --------------------------- a0 * b0, a0 * b1, a0 * b2 end --------------------------*/ \
/* --------------------------- a1 * b0, a1 * b1, a1 * b2 start ------------------------*/ \
/*52*/  "movq "PPS(1,n)"(%%rsi),%%rax             \n" /* a1 into rax */                   \
/*53*/  "movq %%rax            ,%%rcx             \n" /* save a0-rcx faster than stack */ \
/*54*/  "mulq (%%rbx)                             \n" /* a1 * b0 */                       \
/*55*/  "addq %%rax            ,%%r9              \n" /* l46 + a1b0lo */                  \
/*56*/  "adcq %%rdx            ,%%r10             \n" /* l47 + a1b0hi + c */              \
/*57*/  "adcq $0               ,%%r11             \n" /* possible carry, for 192 one adcq 0, 256 two adcq, 320 tree adcq .... */                \
/*58*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a1) from the stack */ \
/*59*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a1*b2 */                         \
/*60*/  "addq %%rax            ,%%r11             \n" /* l57 + a1b2lo + c */              \
/*61*/  "adcq %%rdx            ,%%r12             \n" /* a1b2hi + c */                    \
/*62*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a1) from the stack */ \
/*63*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a1*b1 */                         \
/*64*/  "addq %%rax            ,%%r10             \n" /* a1b2lo to r9 */                  \
/*65*/  "adcq %%rdx            ,%%r11             \n" /* a1b2hi + c  */                   \
/*66*/  "adcq $0               ,%%r12             \n" /* r12 + c  */                      \
/* --------------------------- a1 * b0, a1 * b1, a1 * b2 end --------------------------*/ \
/* --------------------------- a2 * b0, a2 * b1, a2 * b2 start ------------------------*/ \
/*67*/  "movq "PPS(2,n)"(%%rsi),%%rax             \n" /* a2 to rax */                     \
/*68*/  "movq %%rax            ,%%rcx             \n" /* copy into the stack */           \
/*69*/  "mulq (%%rbx)                             \n" /* a2*b0 */                         \
/*70*/  "addq %%rax            ,%%r10             \n" /* l64 + a2b0lo */                  \
/*71*/  "adcq %%rdx            ,%%r11             \n" /* l65 + a2b0hi + c */              \
/*72*/  "adcq $0               ,%%r12             \n" /* possible carry */                \
/*73*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a2) */                \
/*74*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a2*b2 */                         \
/*75*/  "addq %%rax            ,%%r12             \n" /* a2b2lo + l31 + c*/               \
/*76*/  "adcq %%rdx            ,%%r13             \n" /* a2b2hi + l03 + c*/               \
/*77*/  "movq %%rcx            ,%%rax             \n" /* reload a2 */                     \
/*78*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a2*b1 */                         \
/*79*/  "addq %%rax            ,%%r11             \n" /* a2b1lo + l36 */                  \
/*80*/  "adcq %%rdx            ,%%r12             \n" /* a2b1hi + l39 */                  \
/*81*/  "adcq $0               ,%%r13             \n" /* r12 + c */                       \
/* --------------------------- a2 * b0, a2 * b1, a2 * b2 end --------------------------*/ \
/* ---------------------------           sign                --------------------------*/ \
/*81*/  "xorq %%r14            ,%%r15             \n"                                     \
/*82*/  "cmpq $0               ,%%r15             \n" /* r15 = 1 we negate */             \
/*83*/  "je _IsNegativeResult_384_192_            \n" /* not equal ZF = 0, negate*/       \
/*84*/  "notq %%r8                                \n" /* start2ComplementMethod negate */ \
/*85*/  "notq %%r9                                \n" /* 2CM negate */                    \
/*86*/  "notq %%r10                               \n" /* 2CM negate */                    \
/*87*/  "notq %%r11                               \n" /* 2CM negate */                    \
/*88*/  "notq %%r12                               \n" /* 2CM negate */                    \
/*89*/  "notq %%r13                               \n" /* 2CM negate */                    \
/*90*/  "addq $0x1             ,%%r8              \n" /* 2CM add 1 */                     \
/*91*/  "adcq $0x0             ,%%r9              \n" /* 2CM propagate CB */              \
/*92*/  "adcq $0x0             ,%%r10             \n" /* 2CM propagate CB */              \
/*93*/  "adcq $0x0             ,%%r11             \n" /* 2CM propagate CB */              \
/*94*/  "adcq $0x0             ,%%r12             \n" /* 2CM propagate CB */              \
/*95*/  "adcq $0x0             ,%%r13             \n" /* 2CM propagate CB */              \
/*96*/  "_IsNegativeResult_384_192_ :             \n" /* end if negative result */        \
/*97*/  "movq %%r8             ,(%%rdi)           \n" /* r8 -> c1 */                      \
/*98*/  "movq %%r9             ,"PPS(1,n)"(%%rdi) \n" /* r9 -> c1 */                      \
/*99*/  "movq %%r10            ,"PPS(2,n)"(%%rdi) \n" /* r10 -> c2 */                     \
/*100*/ "movq %%r11            ,"PPS(3,n)"(%%rdi) \n" /* r11 -> c3 */                     \
/*101*/ "movq %%r12            ,"PPS(4,n)"(%%rdi) \n" /* r12 -> c4 */                     \
/*102*/ "movq %%r13            ,"PPS(5,n)"(%%rdi) \n" /* r13 -> c5 */                     \
/*103*/ "addq $0x30            ,%%rsp             \n" /* destroy stack frame */           \
   : : : "rax","rdx","rcx","rbx","r8","r9","r10","r11","r12","r13","r14","r15","memory"   \
   ); \
}
// Vli (384) += VLI (192) * VLI (192) 
// c += a*b and than use  operand = sign_bit*(0xFFFFFF) & two_complement | (!sign_bit)*(0xFFFFFFFF) & original_number
#define HELPER_ASM_MUL_ADD384_192_192(n) \
void muladd384_192_192(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */, unsigned long int const* z/* %%rdx -> rbx */){ \
    asm( \
/*-01*/ "subq $0x30            ,%%rsp             \n" /* create stack frame */            \
/*00*/  "movq %%rdx            ,%%rbx             \n" /* rdx uses by mul             */   \
/*01*/  "xorq %%r11            ,%%r11             \n" /* r10 = 0 due to carry effect */   \
/*02*/  "xorq %%r12            ,%%r12             \n" /* r11 = 0 due to carry effect */   \
/*03*/  "xorq %%r13            ,%%r13             \n" /* r12 = 0 due to carry effect */   \
/*03*/  "xorq %%r14            ,%%r14             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/*03*/  "xorq %%r15            ,%%r15             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/* negate a if negative and store into stack, in reverse order due to universal access */ \
/*04*/  "movq "PPS(2,n)"(%%rsi),%%rax             \n" /* load a3 into r8, for the sign */ \
/*05*/  "cmpq $0               ,%%rax             \n" /* test a is negative(sign on a3)*/ \
/*06*/  "jns _MulAddNegativea_384_192_            \n" /* number is negative, we negate */ \
/*07*/  "movq (%%rsi)          ,%%r8              \n" /* load a0 */                       \
/*08*/  "movq "PPS(1,n)"(%%rsi),%%r9              \n" /* load a1 */                       \
/*09*/  "notq %%r8                                \n" /* C2M, ~a0 */                      \
/*10*/  "notq %%r9                                \n" /* C2M, ~a1 */                      \
/*11*/  "notq %%rax                               \n" /* C2M, ~a2 */                      \
/*12*/  "addq $0x1             ,%%r8              \n" /* C2M, ~a0+1 */                    \
/*13*/  "adcq $0x0             ,%%r9              \n" /* C2M, ~a1+CB */                   \
/*14*/  "adcq $0x0             ,%%rax             \n" /* C2M, ~a2+CB */                   \
/*15*/  "movq %%r8             ,-0x18(%%rsp)      \n" /* a0 into the stack -24 rsp */     \
/*16*/  "movq %%r9             ,-0x10(%%rsp)      \n" /* a1 into the stack -16 rsp */     \
/*17*/  "movq %%rax            ,-0x08(%%rsp)      \n" /* a2 into the stack -8 rsp */      \
/*18*/  "leaq  -0x18(%%rsp)    ,%%rsi             \n" /* rsi points to stack a0 > 0 */    \
/*19*/  "movq $1               ,%%r14             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/*20*/  "_MulAddNegativea_384_192_ :              \n" /* end if structure */              \
/* negate a if negative and store into stack, in reverse order due to universal access */ \
/*21*/  "movq "PPS(2,n)"(%%rbx),%%rax             \n" /* load a3 into r8, for the sign */ \
/*22*/  "cmpq $0               ,%%rax             \n" /* test a is negative(sign on a3)*/ \
/*23*/  "jns _MulAddNegativeb_384_192_            \n" /* number is negative, we negate */ \
/*24*/  "movq (%%rbx)          ,%%r8              \n" /* load b0 */                       \
/*25*/  "movq "PPS(1,n)"(%%rbx),%%r9              \n" /* load b1 */                       \
/*26*/  "notq %%r8                                \n" /* C2M, ~b0 */                      \
/*27*/  "notq %%r9                                \n" /* C2M, ~b1 */                      \
/*28*/  "notq %%rax                               \n" /* C2M, ~b2 */                      \
/*29*/  "addq $0x1             ,%%r8              \n" /* C2M, ~b0+1 */                    \
/*30*/  "adcq $0x0             ,%%r9              \n" /* C2M, ~b1+CB */                   \
/*31*/  "adcq $0x0             ,%%rax             \n" /* C2M, ~b2+CB */                   \
/*32*/  "movq %%r8             ,-0x30(%%rsp)      \n" /* b0 into the stack -48 rsp */     \
/*33*/  "movq %%r9             ,-0x28(%%rsp)      \n" /* b1 into the stack -40 rsp */     \
/*34*/  "movq %%rax            ,-0x20(%%rsp)      \n" /* b2 into the stack -32 rsp */     \
/*35*/  "leaq  -0x30(%%rsp)    ,%%rbx             \n" /* rsi points to stack b0 > 0 */    \
/*36*/  "movq $1               ,%%r15             \n" /* r13 = 0 it is the sign 0+ 1-*/   \
/*37*/  "_MulAddNegativeb_384_192_ :              \n" /* end if structure */              \
/*38*/  "xorq %%r10            ,%%r10             \n" /* r9 = 0  due to carry effect */   \
/* --------------------------- a0 * b0, a0 * b1, a0 * b2 start ------------------------*/ \
/*39*/  "movq (%%rsi)          ,%%rax             \n" /* a0 into rax */                   \
/*40*/  "movq %%rax            ,%%rcx             \n" /* save a0-rcx faster than stack */ \
/*41*/  "mulq (%%rbx)                             \n" /* lo rax, hi rdx   a0*b0 */        \
/*42*/  "movq %%rax            ,%%r8              \n" /* only one term, write into c0 */  \
/*43*/  "movq %%rdx            ,%%r9              \n" /* a0b0hi into r8 */                \
/*44*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a0) from the stack */ \
/*45*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a0 * b1 */                       \
/*46*/  "addq %%rax            ,%%r9              \n" /* add a0b0hi + a0b1lo */           \
/*47*/  "adcq %%rdx            ,%%r10             \n" /* save the a0b1hi into r9 */       \
/*48*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a0) from the stack */ \
/*49*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a0 * b2 */                       \
/*50*/  "addq %%rax            ,%%r10             \n" /* add l11 + a0b2lo + c */          \
/*51*/  "adcq %%rdx            ,%%r11             \n" /* add l01 + a0b2hi + c, end */     \
/* --------------------------- a0 * b0, a0 * b1, a0 * b2 end --------------------------*/ \
/* --------------------------- a1 * b0, a1 * b1, a1 * b2 start ------------------------*/ \
/*52*/  "movq "PPS(1,n)"(%%rsi),%%rax             \n" /* a1 into rax */                   \
/*53*/  "movq %%rax            ,%%rcx             \n" /* save a0-rcx faster than stack */ \
/*54*/  "mulq (%%rbx)                             \n" /* a1 * b0 */                       \
/*55*/  "addq %%rax            ,%%r9              \n" /* l46 + a1b0lo */                  \
/*56*/  "adcq %%rdx            ,%%r10             \n" /* l50 + a1b0hi + c */              \
/*57*/  "adcq $0               ,%%r11             \n" /* possible carry */                \
/*58*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a1) from the stack */ \
/*59*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a1*b2 */                         \
/*60*/  "addq %%rax            ,%%r11             \n" /* l57 + a1b2lo + c */              \
/*61*/  "adcq %%rdx            ,%%r12             \n" /* a1b2hi + c */                    \
/*62*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a1) from the stack */ \
/*63*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a1*b1 */                         \
/*64*/  "addq %%rax            ,%%r10             \n" /* a1b2lo to r9 */                  \
/*65*/  "adcq %%rdx            ,%%r11             \n" /* a1b2hi + c  */                   \
/*66*/  "adcq $0               ,%%r12             \n" /* r11 + c  */                      \
/* --------------------------- a1 * b0, a1 * b1, a1 * b2 end --------------------------*/ \
/* --------------------------- a2 * b0, a2 * b1, a2 * b2 start ------------------------*/ \
/*67*/  "movq "PPS(2,n)"(%%rsi),%%rax             \n" /* a2 to rax */                     \
/*68*/  "movq %%rax            ,%%rcx             \n" /* copy into the stack */           \
/*69*/  "mulq (%%rbx)                             \n" /* a2*b0 */                         \
/*70*/  "addq %%rax            ,%%r10             \n" /* l64 + a2b0lo */                  \
/*71*/  "adcq %%rdx            ,%%r11             \n" /* l65 + a2b0hi + c */              \
/*72*/  "adcq $0               ,%%r12             \n" /* possible carry */                \
/*73*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a2) */                \
/*74*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a2*b2 */                         \
/*75*/  "addq %%rax            ,%%r12             \n" /* a2b2lo + l72 + c*/               \
/*76*/  "adcq %%rdx            ,%%r13             \n" /* a2b2hi + l03 + c*/               \
/*77*/  "movq %%rcx            ,%%rax             \n" /* reload a2 */                     \
/*78*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a2*b1 */                         \
/*79*/  "addq %%rax            ,%%r11             \n" /* a2b1lo + l36 */                  \
/*80*/  "adcq %%rdx            ,%%r12             \n" /* a2b1hi + l39 */                  \
/*81*/  "adcq $0               ,%%r13             \n" /* r12 + c */                       \
/* --------------------------- a2 * b0, a2 * b1, a2 * b2 end --------------------------*/ \
/* ---------------------------           sign                --------------------------*/ \
/*81*/  "xorq %%r14            ,%%r15             \n"                                     \
/*82*/  "cmpq $0               ,%%r15             \n" /* r15 = 1 we negate */             \
/*83*/  "je _MulAddIsNegativeResult_384_192_      \n" /* not equal ZF = 0, negate*/       \
/*84*/  "notq %%r8                                \n" /* start2ComplementMethod negate */ \
/*85*/  "notq %%r9                                \n" /* 2CM negate */                    \
/*86*/  "notq %%r10                               \n" /* 2CM negate */                    \
/*87*/  "notq %%r11                               \n" /* 2CM negate */                    \
/*88*/  "notq %%r12                               \n" /* 2CM negate */                    \
/*89*/  "notq %%r13                               \n" /* 2CM negate */                    \
/*90*/  "addq $0x1             ,%%r8              \n" /* 2CM add 1 */                     \
/*91*/  "adcq $0x0             ,%%r9              \n" /* 2CM propagate CB */              \
/*92*/  "adcq $0x0             ,%%r10             \n" /* 2CM propagate CB */              \
/*93*/  "adcq $0x0             ,%%r11             \n" /* 2CM propagate CB */              \
/*94*/  "adcq $0x0             ,%%r12             \n" /* 2CM propagate CB */              \
/*95*/  "adcq $0x0             ,%%r13             \n" /* 2CM propagate CB */              \
/*96*/  "_MulAddIsNegativeResult_384_192_ :       \n" /* end if negative result */        \
/* --------------------------- +=   final addition start     --------------------------*/ \
/*97*/  "addq (%%rdi)          , %%r8             \n" /* add a0+b0 */                     \
/*98*/  "adcq "PPS(1,n)"(%%rdi), %%r9             \n" /* add a1+b1+c */                   \
/*99*/  "adcq "PPS(2,n)"(%%rdi), %%r10            \n" /* add a2+b2+c */                   \
/*100*/ "adcq "PPS(3,n)"(%%rdi), %%r11            \n" /* add a3+b3+c */                   \
/*101*/ "adcq "PPS(4,n)"(%%rdi), %%r12            \n" /* add a4+b4+c */                   \
/*102*/ "adcq "PPS(5,n)"(%%rdi), %%r13            \n" /* add a5+b5+c */                   \
/* --------------------------- +=   end   addition start     --------------------------*/ \
/* ---------------------------  write the result             --------------------------*/ \
/*103*/ "movq %%r8             ,(%%rdi)           \n" /* r8 -> c0 */                      \
/*104*/ "movq %%r9             ,"PPS(1,n)"(%%rdi) \n" /* r8 -> c1 */                      \
/*105*/ "movq %%r10            ,"PPS(2,n)"(%%rdi) \n" /* r9 -> c2 */                      \
/*106*/ "movq %%r11            ,"PPS(3,n)"(%%rdi) \n" /* r10 -> c3 */                     \
/*107*/ "movq %%r12            ,"PPS(4,n)"(%%rdi) \n" /* r11 -> c4 */                     \
/*108*/ "movq %%r13            ,"PPS(5,n)"(%%rdi) \n" /* r12 -> c5 */                     \
/*109*/ "addq $0x30            ,%%rsp             \n" /* destroy stack frame */           \
    : : : "rax","rdx","rbx","r8","r9","r10","r11","r12","r13","r14","r15","memory"        \
    ); \
};

    // c - 1 is the order = AoS, Order*Order = SoA
    // c - generate assembly function 
    HELPER_ASM_MUL384_192_192(1)
    HELPER_ASM_MUL_ADD384_192_192(1)
        
    } //namespase detail
} //namespace vli
