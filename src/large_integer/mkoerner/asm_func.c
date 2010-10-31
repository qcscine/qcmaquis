#include "asm_func.h"

//
// Copyright (c) 2002-2005 Mathias Koerner
//
// Permission is hereby granted, free of charge, to any person obtaining a 
// copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, 
// and/or sell copies of the Software, and to permit persons to whom the 
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
// DEALINGS IN THE SOFTWARE.
//

#ifdef __cplusplus
extern "C" {
#endif

void add_num_num(unsigned long* w, int n, const unsigned long* x) {
  asm("clc
%=:
  movl -4(%%ebx,%%ecx,4),%%edx
  adcl %%edx,-4(%%eax,%%ecx,4)
  decl %%ecx
  jnz %=b"
    : : "c"(n), "a"(w), "b"(x) : "edx", "memory", "cc");
}

void add_num_ui32(unsigned long* w, int n, unsigned long x) {
  asm("addl %2,(%%eax,%%ecx,4)
%=:
  adcl  $0,-4(%%eax,%%ecx,4)
  decl %%ecx
  jnz %=b"
    : : "c"(n-1), "a"(w), "b"(x) : "memory", "cc");
}

void add_num_ui64(unsigned long* w, int n, const unsigned long long x) {
  union {
    unsigned long long i64;
    unsigned long i32[2];
  } tx;
  tx.i64 = x;
  asm("addl %2,4(%%eax,%%ecx,4)
  adcl %3,(%%eax,%%ecx,4)
%=:
  adcl  $0,-4(%%eax,%%ecx,4)
  decl %%ecx # decrement counter
  jnz %=b"
    : : "c"(n-2), "a"(w), "b"(tx.i32[0]), "d"(tx.i32[1]) : "memory", "cc");
}

void sub_num_num(unsigned long* w, int n, const unsigned long* x) {
  asm("clc
%=:
  movl -4(%%ebx,%%ecx,4),%%edx
  sbbl %%edx,-4(%%eax,%%ecx,4)
  decl %%ecx
  jnz %=b"
    : : "c"(n), "a"(w), "b"(x) : "edx", "memory", "cc");
}

void sub_num_ui32(unsigned long* w, int n, unsigned long x) {
  asm("subl %2,(%%eax,%%ecx,4)
%=:
  sbbl  $0,-4(%%eax,%%ecx,4)
  decl %%ecx
  jnz %=b"
    : : "c"(n-1), "a"(w), "b"(x) : "memory", "cc");
}

void sub_num_ui64(unsigned long* w, int n, const unsigned long long x) {
  union {
    unsigned long long i64;
    unsigned long i32[2];
  } tx;
  tx.i64 = x;
  asm("subl %2,4(%%eax,%%ecx,4)
  sbbl %3,(%%eax,%%ecx,4)
%=:
  sbbl  $0,-4(%%eax,%%ecx,4)
  decl %%ecx # decrement counter
  jnz %=b"
    : : "c"(n-2), "a"(w), "b"(tx.i32[0]), "d"(tx.i32[1]) : "memory", "cc");
}

#ifdef __cplusplus
}
#endif
