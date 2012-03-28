//
//  kernels_cpu.cpp
//  VLI_ASM
//
//  Created by Tim Ewart on 22.03.12.
//  Copyright (c) 2012 __Université de Genève__. All rights reserved.
//
#include "kernels_cpu_asm.h"

#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/stringize.hpp>

namespace vli{
    namespace detail{
        
// c- this assembly code support both data layout SoA and AoS
// clear the syntax
#define PPS(m,n) BOOST_PP_STRINGIZE( BOOST_PP_MUL(BOOST_PP_MUL(m,n),8)) // m*n*8, 8 because long int
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

// VLI += long
#define HELPER_ASM_ADD192_64(n) \
void add192_64(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */){ \
    asm( \
        "movq (%%rdi)          , %%r8              \n" /* load a0 */     \
        "movq "PPS(1,n)"(%%rdi), %%r9              \n" /* load a1 */     \
        "movq "PPS(2,n)"(%%rdi), %%rcx             \n" /* load a2 */     \
        "addq (%%rsi)          , %%r8              \n" /* add a0+b0 */   \
        "adcq $0x0             , %%r9              \n" /* add a1+0+c */  \
        "adcq $0x0             , %%rcx             \n" /* add a2+0+c */  \
        "movq %%r8             , (%%rdi)           \n" /* save a0 */     \
        "movq %%r9             , "PPS(1,n)"(%%rdi) \n" /* save a2 */     \
        "movq %%rcx            , "PPS(2,n)"(%%rdi) \n" /* save a2 */     \
        : : :"rcx","r8","r9","memory"                        \
    ); \
};

// VLI += VLI
#define HELPER_ASM_ADD192_192(n) \
void add192_192(unsigned long int* /* %%rdi */, unsigned long int const* /* %%rsi */){ \
    asm( \
        "movq (%%rdi)          , %%r8              \n" /* load a0 */     \
        "movq "PPS(1,n)"(%%rdi), %%r9              \n" /* load a1 */     \
        "movq "PPS(2,n)"(%%rdi), %%rcx             \n" /* load a2 */     \
        "addq (%%rsi)          , %%r8              \n" /* add a0+b0 */   \
        "adcq "PPS(1,n)"(%%rsi), %%r9              \n" /* add a1+b1+c */ \
        "adcq "PPS(2,n)"(%%rsi), %%rcx             \n" /* add a2+b2+c */ \
        "movq %%r8             , (%%rdi)           \n" /* save a0 */     \
        "movq %%r9             , "PPS(1,n)"(%%rdi) \n" /* save a1 */     \
        "movq %%rcx            , "PPS(2,n)"(%%rdi) \n" /* save a2 */     \
        : : :"rcx","r8","r9","memory"                        \
    ); \
};
        
// 2VLI += 2VLI
#define HELPER_ASM_ADD384_384(n) \
void add384_384(unsigned long int* /* %%rdi */, unsigned long int const* /* %%rsi */){ \
    asm( \
        "movq (%%rdi)          , %%r8              \n" /* load a0 */     \
        "movq "PPS(1,n)"(%%rdi), %%r9              \n" /* load a1 */     \
        "movq "PPS(2,n)"(%%rdi), %%r10             \n" /* load a2 */     \
        "movq "PPS(3,n)"(%%rdi), %%r11             \n" /* load a3 */     \
        "movq "PPS(4,n)"(%%rdi), %%r12             \n" /* load a4 */     \
        "movq "PPS(5,n)"(%%rdi), %%rcx             \n" /* load a5 */     \
        "addq (%%rsi)          , %%r8              \n" /* add a0+b0 */   \
        "adcq "PPS(1,n)"(%%rsi), %%r9              \n" /* add a1+b1+c */ \
        "adcq "PPS(2,n)"(%%rsi), %%r10             \n" /* add a2+b2+c */ \
        "adcq "PPS(3,n)"(%%rsi), %%r11             \n" /* add a3+b3+c */ \
        "adcq "PPS(4,n)"(%%rsi), %%r12             \n" /* add a4+b4+c */ \
        "adcq "PPS(5,n)"(%%rsi), %%rcx             \n" /* add a5+b5+c */ \
        "movq %%r8             , (%%rdi)           \n" /* save a0 */     \
        "movq %%r9             , "PPS(1,n)"(%%rdi) \n" /* save a1 */     \
        "movq %%r10            , "PPS(2,n)"(%%rdi) \n" /* save a2 */     \
        "movq %%r11            , "PPS(3,n)"(%%rdi) \n" /* save a3 */     \
        "movq %%r12            , "PPS(4,n)"(%%rdi) \n" /* save a4 */     \
        "movq %%rcx            , "PPS(5,n)"(%%rdi) \n" /* save a5 */     \
        : : :"rcx","r8","r9","r10","r11","r12","memory"                  \
    ); \
};

// 2VLI += long
#define HELPER_ASM_ADD384_64(n) \
void add384_64(unsigned long int* /* %%rdi */, unsigned long int const* /* %%rsi */){ \
    asm( \
        "movq (%%rdi)          , %%r8              \n" /* load a0 */     \
        "movq "PPS(1,n)"(%%rdi), %%r9              \n" /* load a1 */     \
        "movq "PPS(2,n)"(%%rdi), %%r10             \n" /* load a2 */     \
        "movq "PPS(3,n)"(%%rdi), %%r11             \n" /* load a3 */     \
        "movq "PPS(4,n)"(%%rdi), %%r12             \n" /* load a4 */     \
        "movq "PPS(5,n)"(%%rdi), %%rcx             \n" /* load a5 */     \
        "addq (%%rsi)          , %%r8              \n" /* add a0+b0 */   \
        "adcq $0x0             , %%r9              \n" /* add a1+b1+c */ \
        "adcq $0x0             , %%r10             \n" /* add a2+0+c */  \
        "adcq $0x0             , %%r11             \n" /* add a3+0+c */  \
        "adcq $0x0             , %%r12             \n" /* add a4+0+c */  \
        "adcq $0x0             , %%rcx             \n" /* add a5+0+c */  \
        "movq %%r8             , (%%rdi)           \n" /* save a0 */     \
        "movq %%r9             , "PPS(1,n)"(%%rdi) \n" /* save a1 */     \
        "movq %%r10            , "PPS(2,n)"(%%rdi) \n" /* save a2 */     \
        "movq %%r11            , "PPS(3,n)"(%%rdi) \n" /* save a3 */     \
        "movq %%r12            , "PPS(4,n)"(%%rdi) \n" /* save a4 */     \
        "movq %%rcx            , "PPS(5,n)"(%%rdi) \n" /* save a5 */     \
        : : :"rcx","r8","r9","r10","r11","r12","memory"                  \
    ); \
};

//VLI -= long
#define HELPER_ASM_SUB192_64(n) \
void sub192_64(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */){ \
    asm( \
        "movq (%%rdi)          , %%r8              \n" /* load a0 */     \
        "movq "PPS(1,n)"(%%rdi), %%r9              \n" /* load a0 */     \
        "movq "PPS(2,n)"(%%rdi), %%rcx             \n" /* load a0 */     \
        "subq (%%rsi)          , %%r8              \n" /* sub a0-b0 */   \
        "sbbq $0x0             , %%r9              \n" /* sub a1-0-b */  \
        "sbbq $0x0             , %%rcx             \n" /* sub a2-0-b */  \
        "movq %%r8             , (%%rdi)           \n" /* save a0 */     \
        "movq %%r9             , "PPS(1,n)"(%%rdi) \n" /* save a1 */     \
        "movq %%rcx            , "PPS(2,n)"(%%rdi) \n" /* save a2 */     \
        : : :"rcx","r8","r9","memory"                        \
    ); \
};

// VLI -= VLI
#define HELPER_ASM_SUB192_192(n) \
void sub192_192(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */){ \
    asm( \
        "movq (%%rdi)          , %%r8              \n" /* load a0 */     \
        "movq "PPS(1,n)"(%%rdi), %%r9              \n" /* load a1 */     \
        "movq "PPS(2,n)"(%%rdi), %%rcx             \n" /* load a2 */     \
        "subq (%%rsi)          , %%r8              \n" /* sub a0-b0 */   \
        "sbbq "PPS(1,n)"(%%rsi), %%r9              \n" /* sub a1-b1-b */ \
        "sbbq "PPS(2,n)"(%%rsi), %%rcx             \n" /* sub a2-b2-b */ \
        "movq %%r8             , (%%rdi)           \n" /* save a0 */     \
        "movq %%r9             , "PPS(1,n)"(%%rdi) \n" /* save a1 */     \
        "movq %%rcx            , "PPS(2,n)"(%%rdi) \n" /* save a2 */     \
        : : :"rcx","r8","r9","memory"                        \
    ); \
};

// 2*VLI -= long 
#define HELPER_ASM_SUB384_64(n) \
void sub384_64(unsigned long int* /* %%rdi */, unsigned long int const* /* %%rsi */){ \
    asm( \
        "movq (%%rdi)          , %%r8              \n" /* load a0 */     \
        "movq "PPS(1,n)"(%%rdi), %%r9              \n" /* load a1 */     \
        "movq "PPS(2,n)"(%%rdi), %%r10             \n" /* load a2 */     \
        "movq "PPS(3,n)"(%%rdi), %%r11             \n" /* load a3 */     \
        "movq "PPS(4,n)"(%%rdi), %%r12             \n" /* load a4 */     \
        "movq "PPS(5,n)"(%%rdi), %%rcx             \n" /* load a5 */     \
        "subq (%%rsi)          , %%r8              \n" /* add a0-b0 */   \
        "sbbq $0x0             , %%r9              \n" /* add a1-b1-c */ \
        "sbbq $0x0             , %%r10             \n" /* add a2-c */    \
        "sbbq $0x0             , %%r11             \n" /* add a3-c */    \
        "sbbq $0x0             , %%r12             \n" /* add a4-c */    \
        "sbbq $0x0             , %%rcx             \n" /* add a5-c */    \
        "movq %%r8             , (%%rdi)           \n" /* save a0 */     \
        "movq %%r9             , "PPS(1,n)"(%%rdi) \n" /* save a1 */     \
        "movq %%r10            , "PPS(2,n)"(%%rdi) \n" /* save a2 */     \
        "movq %%r11            , "PPS(3,n)"(%%rdi) \n" /* save a3 */     \
        "movq %%r12            , "PPS(4,n)"(%%rdi) \n" /* save a4 */     \
        "movq %%rcx            , "PPS(5,n)"(%%rdi) \n" /* save a5 */     \
        : : :"rcx","r8","r9","r10","r11","r12","memory"                  \
    ); \
};

// 2*VLI -= 2*VLI
#define HELPER_ASM_SUB384_384(n) \
void sub384_384(unsigned long int* /* %%rdi */, unsigned long int const* /* %%rsi */){ \
    asm( \
        "movq (%%rdi)          , %%r8              \n" /* load a0 */     \
        "movq "PPS(1,n)"(%%rdi), %%r9              \n" /* load a1 */     \
        "movq "PPS(2,n)"(%%rdi), %%r10             \n" /* load a2 */     \
        "movq "PPS(3,n)"(%%rdi), %%r11             \n" /* load a3 */     \
        "movq "PPS(4,n)"(%%rdi), %%r12             \n" /* load a4 */     \
        "movq "PPS(5,n)"(%%rdi), %%rcx             \n" /* load a5 */     \
        "subq (%%rsi)          , %%r8              \n" /* add a0-b0 */   \
        "sbbq "PPS(1,n)"(%%rsi), %%r9              \n" /* add a1-b1-c */ \
        "sbbq "PPS(2,n)"(%%rsi), %%r10             \n" /* add a2-b2-c */ \
        "sbbq "PPS(3,n)"(%%rsi), %%r11             \n" /* add a3-b3-c */ \
        "sbbq "PPS(4,n)"(%%rsi), %%r12             \n" /* add a4-b4-c */ \
        "sbbq "PPS(5,n)"(%%rsi), %%rcx             \n" /* add a5-b5-c */ \
        "movq %%r8             , (%%rdi)           \n" /* save a0 */     \
        "movq %%r9             , "PPS(1,n)"(%%rdi) \n" /* save a1 */     \
        "movq %%r10            , "PPS(2,n)"(%%rdi) \n" /* save a2 */     \
        "movq %%r11            , "PPS(3,n)"(%%rdi) \n" /* save a3 */     \
        "movq %%r12            , "PPS(4,n)"(%%rdi) \n" /* save a4 */     \
        "movq %%rcx            , "PPS(5,n)"(%%rdi) \n" /* save a5 */     \
        : : :"rcx","r8","r9","r10","r11","r12","memory"                  \
    ); \
};

// Vli *= long
#define HELPER_ASM_MUL192_64(n) \
void mul192_64(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */){ \
    asm( \
        "movq (%%rsi)          ,%%rax             \n" /* a0 into rax */                   \
        "movq %%rbp            ,-0x8(%%rsp)       \n" /* alloc the stack */               \
        "xorq %%rcx            ,%%rcx             \n" /* rcx to 0 */                      \
        "cmpq %%rax            ,%%rcx             \n" /* rax is negative ? */             \
        "js   _IsNegative192_64                   \n" /* yes it is, SF = 1 */             \
        "negq %%rax                               \n" /* negate the number */             \
        "movq $0x1             ,%%rcx             \n" /* keep trace for final sign */     \
        "_IsNegative192_64 :                      \n" /* end if structure */              \
        "movq %%rax            ,%%rbp             \n" /* keep a copy of rax/a0 inside */  \
        "mulq (%%rdi)                             \n" /* lo rax, hi rdx   a0*b0 */        \
        "movq %%rax            ,%%r8              \n" /* only one term, write into c0 */  \
        "movq %%rdx            ,%%r9              \n" /* hia0b0 into rcx */               \
        "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
        "mulq "PPS(1,n)"(%%rdi)                   \n" /* a0 * b1 */                       \
        "addq %%rax            ,%%r9              \n" /* add hia0b0 + loa0b1 */           \
        "movq %%rdx            ,%%r10             \n" /* save the hi into rcx */          \
        "adcq $0               ,%%r10             \n" /* perhaps carry */                 \
        "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
        "imulq "PPS(2,n)"(%%rdi)                  \n" /* a0 * b2, we skip the the hi */   \
        "addq %%rax            ,%%r10             \n" /* add hia0b1 + loa0b2 */           \
        "cmpq $0               ,%%rcx             \n" /* rcx = 1 we negate */             \
        "je _IsNegativeResult192_64               \n" /* not equal ZF = 0, negate*/       \
        "notq %%r8                                \n" /* start C2M negate */              \
        "notq %%r9                                \n" /* 2ComplementMethod negate */      \
        "notq %%r10                               \n" /* 2CM negate */                    \
        "addq $0x1             ,%%r8              \n" /* 2CM add 1 */                     \
        "adcq $0x0             ,%%r9              \n" /* 2CM propagate CB */              \
        "adcq $0x0             ,%%r10             \n" /* 2CM propagate CB */              \
        "_IsNegativeResult192_64 :                \n" /* end if*/                         \
        "movq %%r8             ,(%%rdi)           \n" /* move into a0 */                  \
        "movq %%r9             ,"PPS(1,n)"(%%rdi) \n" /* move into a1 */                  \
        "movq %%r10            ,"PPS(2,n)"(%%rdi) \n" /* move into a2 */                  \
        "movq -0x8(%%rsp)      ,%%rbp             \n" /* stack clean up */                \
        : : :"rax","rdx","rcx","r8","r9","memory"                 \
    ); \
}

// 2*Vli *= long
#define HELPER_ASM_MUL384_64(n) \
void mul384_64(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */){ \
    asm( \
        "movq (%%rsi)          ,%%rax             \n" /* a0 into rax */                   \
        "movq %%rbp            ,-0x8(%%rsp)       \n" /* alloc the stack */               \
        "xorq %%rcx            ,%%rcx             \n" /* rcx to 0 */                      \
        "cmpq %%rax            ,%%rcx             \n" /* rax is negative ? */             \
        "js   _IsNegative384_64                   \n" /* yes it is, SF = 1 */             \
        "negq %%rax                               \n" /* negate the number */             \
        "movq $0x1             ,%%rcx             \n" /* keep trace for final sign */     \
        "_IsNegative384_64 :                      \n" /* end if structure */              \
        "movq %%rax            ,%%rbp             \n" /* keep a copy of rax/a0 inside */  \
        "mulq (%%rdi)                             \n" /* lo rax, hi rdx   a0*b0 */        \
        "movq %%rax            ,%%r8              \n" /* only one term, write into b0 */  \
        "movq %%rdx            ,%%r9              \n" /* hia0b0 into rcx */               \
        "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
        "mulq "PPS(1,n)"(%%rdi)                   \n" /* a0 * b1 */                       \
        "addq %%rax            ,%%r9              \n" /* add hia0b1 + loa0b1 */           \
        "movq %%rdx            ,%%r10             \n" /* save the hi into rcx */          \
        "adcq $0               ,%%r10             \n" /* perhaps carry */                 \
        "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
        "mulq "PPS(2,n)"(%%rdi)                   \n" /* a0 * b2 */                       \
        "addq %%rax            ,%%r10             \n" /* add hia0b2 + loa0b2 */           \
        "movq %%rdx            ,%%r11             \n" /* save the hi into rcx */          \
        "adcq $0               ,%%r11             \n" /* perhaps carry */                 \
        "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
        "mulq "PPS(3,n)"(%%rdi)                   \n" /* a0 * b3 */                       \
        "addq %%rax            ,%%r11             \n" /* add hia0b3 + loa0b3 */           \
        "movq %%rdx            ,%%r12             \n" /* save the hi into rcx */          \
        "adcq $0               ,%%r12             \n" /* perhaps carry */                 \
        "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
        "mulq "PPS(4,n)"(%%rdi)                   \n" /* a0 * b4 */                       \
        "addq %%rax            ,%%r12             \n" /* add hia0b4 + loa0b4 */           \
        "movq %%rdx            ,%%r13             \n" /* save the hi into rcx */          \
        "adcq $0               ,%%r13             \n" /* perhaps carry */                 \
        "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
        "imulq "PPS(5,n)"(%%rdi)                  \n" /* a0 * b5, we skip the the hi */   \
        "addq %%rax            ,%%r13             \n" /* add hia0b5 + loa0b5 */           \
        "cmpq $0               ,%%rcx             \n" /* rcx = 1 we negate */             \
        "je _IsNegativeResult384_64               \n" /* not equal ZF = 0, negate*/       \
        "notq %%r8                                \n" /* start2ComplementMethod negate */ \
        "notq %%r9                                \n" /* 2CM negate */                    \
        "notq %%r10                               \n" /* 2CM negate */                    \
        "notq %%r11                               \n" /* 2CM negate */                    \
        "notq %%r12                               \n" /* 2CM negate */                    \
        "notq %%r13                               \n" /* 2CM negate */                    \
        "addq $0x1             ,%%r8              \n" /* 2CM add 1 */                     \
        "adcq $0x0             ,%%r9              \n" /* 2CM propagate CB */              \
        "adcq $0x0             ,%%r10             \n" /* 2CM propagate CB */              \
        "adcq $0x0             ,%%r11             \n" /* 2CM propagate CB */              \
        "adcq $0x0             ,%%r12             \n" /* 2CM propagate CB */              \
        "adcq $0x0             ,%%r13             \n" /* 2CM propagate CB */              \
        "_IsNegativeResult384_64 :                \n" /* end if*/                         \
        "movq %%r8             ,(%%rdi)           \n" /* move into a0 */                  \
        "movq %%r9             ,"PPS(1,n)"(%%rdi) \n" /* move into a1 */                  \
        "movq %%r10            ,"PPS(2,n)"(%%rdi) \n" /* move into a2 */                  \
        "movq %%r11            ,"PPS(3,n)"(%%rdi) \n" /* move into a2 */                  \
        "movq %%r12            ,"PPS(4,n)"(%%rdi) \n" /* move into a2 */                  \
        "movq %%r13            ,"PPS(5,n)"(%%rdi) \n" /* move into a2 */                  \
        "movq -0x8(%%rsp)      ,%%rbp             \n" /* stack clean up */                \
        : : :"rax","rdx","rcx","r8","r9","r10","r11","r12","r13","memory"                 \
    ); \
}

// Vli (192) = VLI (192) * VLI (192) 
#define HELPER_ASM_MUL192_192(n) \
void mul192_192(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */){ \
     asm( \
/*01*/  "xorq %%r10            ,%%r10             \n" /* r10 = 0 due to carry effect */   \
/*02*/  "movq (%%rsi)          ,%%rax             \n" /* a0 into rax */                   \
/*03*/  "movq %%rbp            ,-0x08(%%rsp)      \n" /* alloc the stack */               \
/*04*/  "movq %%rax            ,%%rbp             \n" /* keep a stack copy of rax/a0 */   \
/* --------------------------- a0 * b0, a0 * b1, a0 * b2 start ------------------------*/ \
/*05*/  "mulq (%%rdi)                             \n" /* lo rax, hi rdx   a0*b0 */        \
/*06*/  "movq %%rax            ,%%r8              \n" /* only one term, write into c0 */  \
/*07*/  "movq %%rdx            ,%%r9              \n" /* a0b0hi into r8 */                \
/*08*/  "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
/*09*/  "mulq "PPS(1,n)"(%%rdi)                   \n" /* a0 * b1 */                       \
/*10*/  "addq %%rax            ,%%r9              \n" /* add a0b0hi + a0b1lo */           \
/*11*/  "adcq %%rdx            ,%%r10             \n" /* save the a0b1hi into r9 */       \
/*12*/  "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
/*13*/  "imulq "PPS(2,n)"(%%rdi)                  \n" /* a0 * b2 */                       \
/*14*/  "addq %%rax            ,%%r10             \n" /* add l11 + a0b2lo + c */          \
/* --------------------------- a0 * b0, a0 * b1, a0 * b2 end --------------------------*/ \
/* --------------------------- a1 * b0, a1 * b start ----------------------------------*/ \
/*15*/  "movq "PPS(1,n)"(%%rsi),%%rax             \n" /* a1 into rax */                   \
/*16*/  "movq %%rax            ,%%rbp             \n" /* keep a stack copy of rax/a1 */   \
/*17*/  "mulq (%%rdi)                             \n" /* a1 * b0 */                       \
/*18*/  "addq %%rax            ,%%r9              \n" /* l13 + a1b0lo */                  \
/*19*/  "adcq %%rdx            ,%%r10             \n" /* l17 + a1b0hi + c */              \
/*20*/  "movq %%rbp            ,%%rax             \n" /* reload rax(a1) from the stack */ \
/*21*/  "imulq "PPS(1,n)"(%%rdi)                  \n" /* a1*b1 */                         \
/*22*/  "addq %%rax            ,%%r10             \n" /* a1b2lo to r10 */                 \
/* --------------------------- a1 * b0, a1 * b1, a1 * b2 end --------------------------*/ \
/* --------------------------- a2 * b0                   start ------------------------*/ \
/*23*/  "movq "PPS(2,n)"(%%rsi),%%rax             \n" /* a2 to rax */                     \
/*25*/  "imulq (%%rdi)                            \n" /* a2*b0 */                         \
/*26*/  "addq %%rax            ,%%r10             \n" /* l30 + a2b0lo */                  \
/* --------------------------- a2 * b0                   end --------------------------*/ \
/*27*/  "movq %%r8             ,(%%rdi)           \n" /* r8 -> c0 */                      \
/*28*/  "movq %%r9             ,"PPS(1,n)"(%%rdi) \n" /* r9 -> c1 */                      \
/*28*/  "movq %%r10            ,"PPS(2,n)"(%%rdi) \n" /* r10 -> c2 */                     \
/*29*/  "movq -0x08(%%rsp)      ,%%rbp            \n" /* stack clean up */                \
        : : :"rax","rdx","r8","r9","r10","memory" \
    ); \
};

// Vli (384) = VLI (192) * VLI (192) 
#define HELPER_ASM_MUL384_192_192(n) \
void mul384_192_192(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */, unsigned long int const* z/* %%rdx -> rbx */){ \
   asm( \
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
/*15*/  "movq %%r8             ,-24(%%rsp)        \n" /* a0 into the stack -24 rsp */     \
/*16*/  "movq %%r9             ,-16(%%rsp)        \n" /* a1 into the stack -16 rsp */     \
/*17*/  "movq %%rax            ,-8(%%rsp)         \n" /* a2 into the stack -8 rsp */      \
/*18*/  "leaq  -24(%%rsp)      ,%%rsi             \n" /* rsi points to stack a0 > 0 */    \
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
/*32*/  "movq %%r8             ,-48(%%rsp)        \n" /* b0 into the stack -48 rsp */     \
/*33*/  "movq %%r9             ,-40(%%rsp)        \n" /* b1 into the stack -40 rsp */     \
/*34*/  "movq %%rax            ,-32(%%rsp)        \n" /* b2 into the stack -32 rsp */     \
/*35*/  "leaq  -48(%%rsp)      ,%%rbx             \n" /* rsi points to stack b0 > 0 */    \
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
/*55*/  "addq %%rax            ,%%r9              \n" /* l13 + a1b0lo */                  \
/*56*/  "adcq %%rdx            ,%%r10             \n" /* l17 + a1b0hi + c */              \
/*57*/  "adcq $0               ,%%r11             \n" /* possible carry */                \
/*58*/  "movq %%rcx            ,%%rax             \n" /* reload rax(a1) from the stack */ \
/*59*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a1*b2 */                         \
/*60*/  "addq %%rax            ,%%r11             \n" /* l18 + a1b2lo + c */              \
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
/*70*/  "addq %%rax            ,%%r10             \n" /* l30 + a2b0lo */                  \
/*71*/  "adcq %%rdx            ,%%r11             \n" /* l31 + a2b0hi + c */              \
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
/*84*/   "notq %%r8                               \n" /* start2ComplementMethod negate */ \
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
   : : : "rax","rdx","rcx","rbx","r8","r9","r10","r11","r12","r13","r14","r15","memory"   \
   ); \
}
// Vli (384) = VLI (192) * VLI (192) 
// c += a*b
#define HELPER_ASM_MUL_ADD384_192_192(n) \
void muladd384_192_192(unsigned long int* x/* %%rdi */, unsigned long int const* y/* %%rsi */, unsigned long int const* z/* %%rdx -> rbx */){ \
    asm( \
/* ----------- * we alloc into the register for the intermediate result */                \
/*00*/  "movq %%rdx            ,%%rbx             \n" /* rdx uses by mul  */              \
/*01*/  "xorq %%r8             ,%%r8              \n" /* r8 = 0  */                       \
/*02*/  "xorq %%r9             ,%%r9              \n" /* r9 = 0  */                       \
/*03*/  "xorq %%r10            ,%%r10             \n" /* r10 = 0 */                       \
/*04*/  "xorq %%r11            ,%%r11             \n" /* r11 = 0 */                       \
/*05*/  "xorq %%r12            ,%%r12             \n" /* r12 = 0 */                       \
/*06*/  "xorq %%r13            ,%%r13             \n" /* r12 = 0 */                       \
/*07*/  "movq (%%rsi)          ,%%rax             \n" /* a0 into rax */                   \
/*08*/  "movq %%rbp            ,-0x08(%%rsp)      \n" /* alloc the stack */               \
/*09*/  "movq %%rax            ,%%rbp             \n" /* keep a stack copy of rax/a0 */   \
/* --------------------------- a0 * b0, a0 * b1, a0 * b2 start ------------------------*/ \
/*10*/  "mulq (%%rbx)                             \n" /* lo rax, hi rdx   a0*b0 */        \
/*11*/  "movq %%rax            ,%%r8              \n" /* only one term, write into c0 */  \
/*12*/  "movq %%rdx            ,%%r9              \n" /* a0b0hi into r9 */                \
/*13*/  "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
/*14*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a0 * b1 */                       \
/*15*/  "addq %%rax            ,%%r9              \n" /* add a0b0hi + a0b1lo */           \
/*16*/  "adcq %%rdx            ,%%r10             \n" /* save the a0b1hi into r10 */      \
/*17*/  "movq %%rbp            ,%%rax             \n" /* reload rax(a0) from the stack */ \
/*18*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a0 * b2 */                       \
/*19*/  "addq %%rax            ,%%r10             \n" /* add l16 + a0b2lo + c */          \
/*20*/  "adcq %%rdx            ,%%r11             \n" /* add l04 + a0b2hi + c, end */     \
/* --------------------------- a0 * b0, a0 * b1, a0 * b2 end --------------------------*/ \
/* --------------------------- a1 * b0, a1 * b1, a1 * b2 start ------------------------*/ \
/*21*/  "movq "PPS(1,n)"(%%rsi),%%rax             \n" /* a1 into rax */                   \
/*22*/  "movq %%rax            ,%%rbp             \n" /* keep a stack copy of rax/a1 */   \
/*23*/  "mulq (%%rbx)                             \n" /* a1 * b0 */                       \
/*24*/  "addq %%rax            ,%%r9              \n" /* l15 + a1b0lo */                  \
/*25*/  "adcq %%rdx            ,%%r10             \n" /* l16 + a1b0hi + c */              \
/*26*/  "adcq $0               ,%%r11             \n" /* possible carry */                \
/*27*/  "movq %%rbp            ,%%rax             \n" /* reload rax(a1) from the stack */ \
/*28*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a1*b2 */                         \
/*29*/  "addq %%rax            ,%%r11             \n" /* l26 + a1b2lo + c */              \
/*30*/  "adcq %%rdx            ,%%r12             \n" /* a1b2hi + c */                    \
/*31*/  "movq %%rbp            ,%%rax             \n" /* reload rax(a1) from the stack */ \
/*32*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a1*b1 */                         \
/*33*/  "addq %%rax            ,%%r10             \n" /* a1b2lo to r10 */                 \
/*34*/  "adcq %%rdx            ,%%r11             \n" /* a1b2hi + c  */                   \
/*35*/  "adcq $0               ,%%r12             \n" /* r12 + c  */                      \
/* --------------------------- a1 * b0, a1 * b1, a1 * b2 end --------------------------*/ \
/* --------------------------- a2 * b0, a2 * b1, a2 * b2 start ------------------------*/ \
/*36*/  "movq "PPS(2,n)"(%%rsi),%%rax             \n" /* a2 to rax */                     \
/*37*/  "movq %%rax            ,%%rbp             \n" /* copy into the stack */           \
/*38*/  "mulq (%%rbx)                             \n" /* a2*b0 */                         \
/*39*/  "addq %%rax            ,%%r10             \n" /* l33 + a2b0lo */                  \
/*40*/  "adcq %%rdx            ,%%r11             \n" /* l34 + a2b0hi + c */              \
/*41*/  "adcq $0               ,%%r12             \n" /* possible carry */                \
/*42*/  "movq %%rbp            ,%%rax             \n" /* reload rax(a2) */                \
/*43*/  "mulq "PPS(2,n)"(%%rbx)                   \n" /* a2*b2 */                         \
/*44*/  "addq %%rax            ,%%r12             \n" /* a2b2lo + l41 + c*/               \
/*45*/  "adcq %%rdx            ,%%r13             \n" /* a2b2hi + l13 + c*/               \
/*46*/  "movq %%rbp            ,%%rax             \n" /* reload a2 */                     \
/*47*/  "mulq "PPS(1,n)"(%%rbx)                   \n" /* a2*b1 */                         \
/*48*/  "addq %%rax            ,%%r11             \n" /* a2b1lo + l40 */                  \
/*49*/  "adcq %%rdx            ,%%r12             \n" /* a2b1hi + l44 */                  \
/*50*/  "adcq $0               ,%%r13             \n" /* r12 + c */                       \
/* --------------------------- a2 * b0, a2 * b1, a2 * b2 end --------------------------*/ \
/* --------------------------- +=   final addition start     --------------------------*/ \
/*51*/  "addq (%%rdi)          , %%r8             \n" /* add a0+b0 */                     \
/*52*/  "adcq "PPS(1,n)"(%%rdi), %%r9             \n" /* add a1+b1+c */                   \
/*53*/  "adcq "PPS(2,n)"(%%rdi), %%r10            \n" /* add a2+b2+c */                   \
/*54*/  "adcq "PPS(3,n)"(%%rdi), %%r11            \n" /* add a3+b3+c */                   \
/*55*/  "adcq "PPS(4,n)"(%%rdi), %%r12            \n" /* add a4+b4+c */                   \
/*56*/  "adcq "PPS(5,n)"(%%rdi), %%r13            \n" /* add a5+b5+c */                   \
/* --------------------------- +=   end   addition start     --------------------------*/ \
/* ---------------------------  write the result             --------------------------*/ \
/*57*/  "movq %%r8             ,(%%rdi)           \n" /* r8 -> c0 */                      \
/*58*/  "movq %%r9             ,"PPS(1,n)"(%%rdi) \n" /* r8 -> c1 */                      \
/*59*/  "movq %%r10            ,"PPS(2,n)"(%%rdi) \n" /* r9 -> c2 */                      \
/*60*/  "movq %%r11            ,"PPS(3,n)"(%%rdi) \n" /* r10 -> c3 */                     \
/*61*/  "movq %%r12            ,"PPS(4,n)"(%%rdi) \n" /* r11 -> c4 */                     \
/*62*/  "movq %%r13            ,"PPS(5,n)"(%%rdi) \n" /* r12 -> c5 */                     \
/*63*/  "movq -0x08(%%rsp)     ,%%rbp             \n" /* stack clean up */                \
    : : : "rax","rdx","rbx","r8","r9","r10","r11","r12","r13","memory"  \
    ); \
};

    // c - 1 is the order = AoS, Order*Order = SoA
    // c - generate assembly function 
    HELPER_ASM_ADD192_64(1) 
    HELPER_ASM_ADD192_192(1)
    HELPER_ASM_ADD384_64(1)
    HELPER_ASM_ADD384_384(1)
    HELPER_ASM_SUB192_64(1)
    HELPER_ASM_SUB192_192(1)
    HELPER_ASM_SUB384_64(1)
    HELPER_ASM_SUB384_384(1)
    HELPER_ASM_MUL192_64(1)
    HELPER_ASM_MUL192_192(1)
    HELPER_ASM_MUL384_64(1)
    HELPER_ASM_MUL384_192_192(1)
    HELPER_ASM_MUL_ADD384_192_192(1)
        
    } //namespase detail
} //namespace vli
