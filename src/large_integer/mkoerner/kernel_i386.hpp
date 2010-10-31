#ifndef KERNEL_I386_HPP
#define KERNEL_I386_HPP

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

#include <cstddef>
// deque is used for output.
// it should be moved to a different file at some point.
#include <deque>

namespace i386_asm_kernel {

  // typedefs of the basic arithmetic types we
  // are using.
  typedef long int32_t;
  typedef unsigned long uint32_t;
  typedef long long int64_t;
  typedef unsigned long long uint64_t;

  template <typename T> struct closest_type
  { typedef uint32_t type; };
  template <> struct closest_type<int>
  { typedef uint32_t type; };
  template <> struct closest_type<unsigned int>
  { typedef uint32_t type; };
  template <> struct closest_type<long>
  { typedef uint32_t type; };
  template <> struct closest_type<long long>
  { typedef uint64_t type; };

  // a small union that helps us to extract the high
  // and low dword from a quadword.
  // i32[0] is the low dword, i32[1] the high dword.
  union uint64_hilo {
    uint64_t i64;
    uint32_t i32[2];
  };

  template <std::size_t N>
  struct kernel {

    struct num_t {
      uint32_t data[N];
    };

    static inline void init(num_t&);
    static inline void init(num_t&, const uint32_t&);
    static inline void init(num_t&, const uint64_t&);
    static inline void init(num_t&, const num_t&);
    static inline void init(num_t&, const char*);
    template <typename T>
    static inline void init(num_t& t, const T& x)
    { init(t,static_cast<typename closest_type<T>::type>(x)); }

    template <class OSTREAM>
    static void print(const num_t&, OSTREAM&);

    static inline void add(num_t&, const num_t&);
    static inline void add(num_t&, const uint32_t&);
    static inline void add(num_t&, const uint64_t&);
    template <typename T>
    static inline void add(num_t& t, const T& x)
    { add(t,static_cast<typename closest_type<T>::type>(x)); }

    static inline void sub(num_t&, const num_t&);
    static inline void sub(num_t&, const uint32_t&);
    static inline void sub(num_t&, const uint64_t&);
    template <typename T>
    static inline void sub(num_t& t, const T& x)
    { sub(t,static_cast<typename closest_type<T>::type>(x)); }

    static inline void mul(num_t&, const num_t&, const num_t&);
    static inline void mul(num_t&, const num_t&, const uint32_t&);
    static inline void mul(num_t&, const num_t&, const uint64_t&);
    template <typename T>
    static inline void mul(num_t& t, const num_t& t2, const T& x)
    { mul(t,t2,static_cast<typename closest_type<T>::type>(x)); }

    static inline void div(num_t&, uint32_t&, const num_t&, const uint32_t&);

    static inline bool less(const num_t&, const num_t&);
    static inline bool less(const num_t&, const uint32_t&);
    static inline bool less(const num_t&, const uint64_t&);
    template <typename T>
    static inline bool less(const num_t& t, const T& x)
    { return less(t,static_cast<typename closest_type<T>::type>(x)); }

    static inline bool equal(const num_t&, const num_t&);
    static inline bool equal(const num_t&, const uint32_t&);
    static inline bool equal(const num_t&, const uint64_t&);
    template <typename T>
    static inline bool equal(const num_t& t, const T& x)
    { return equal(t,static_cast<typename closest_type<T>::type>(x)); }

    static inline bool greater(const num_t&, const num_t&);
    static inline bool greater(const num_t&, const uint32_t&);
    static inline bool greater(const num_t&, const uint64_t&);
    template <typename T>
    static inline bool greater(const num_t& t, const T& x)
    { return greater(t,static_cast<typename closest_type<T>::type>(x)); }

  };

  inline void to_base10(std::deque<char>&, uint32_t);

  ////////////////////////IMPLEMENTATION/////////////////////////

  // helper function to print numbers in base 10
  void to_base10(std::deque<char>& str, uint32_t l) {
    const int digits_per_ui32 = 9; // 10^9 < 2^32 < 10^10
    for ( int i = 0; i < digits_per_ui32; ++i ) {
      char digit = '0' + ( l % 10 );
      l /= 10;
      str.push_front(digit);
    }
  }

  template <std::size_t N>
  void kernel<N>::init(num_t& x) {
    for ( int i = 0; i < N; ++i )
      x.data[i] = 0UL;
  }

  template <std::size_t N>
  void kernel<N>::init(num_t& x, const uint32_t& y) {
    for ( int i = 0; i < N-1; ++i )
      x.data[i] = 0UL;
    x.data[N-1] = y;
  }

  template <std::size_t N>
  void kernel<N>::init(num_t& x, const uint64_t& y) {
    uint64_hilo tmp;
    tmp.i64 = y;
    for ( int i = 0; i < N-2; ++i )
      x.data[i] = 0U;
    x.data[N-1] = tmp.i32[0];
    x.data[N-2] = tmp.i32[1];
  }

  template <std::size_t N>
  void kernel<N>::init(num_t& x, const num_t& y) {
    for ( int i = 0; i < N; ++i )
      x.data[i] = y.data[i];
  }

  template <std::size_t N>
  void kernel<N>::init(num_t& x, const char* y) {
    init(x);
    uint32_t base = 10UL;
    while ( (*y) != 0UL ) {
      uint32_t digit = (*y) - '0';
      num_t tmp;
      mul(tmp,x,base);
      add(tmp,digit);
      init(x,tmp);
      ++y;
    }
  }

  template <std::size_t N> template <class OSTREAM>
  void kernel<N>::print(const num_t& x, OSTREAM& os) {
    if ( equal(x,0UL) )
      os << 0UL;
    else {
      std::deque<char> str;
      num_t tmp, q;
      init(tmp,x);
      do {
  uint32_t r;
  div(q,r,tmp,1000000000UL);
  init(tmp,q);
  to_base10(str,r);
      } while ( !equal(q,0UL) );
      std::deque<char>::const_iterator i, iend;
      i = str.begin();
      iend = str.end();
      // skip leading 0
      while ( *i == '0' )
  ++i;
      for ( ; i != iend; ++i )
  os << *i;
    }
  }

  template <std::size_t N>
  void kernel<N>::add(num_t& x, const num_t& y) {
    asm("movl %0,%%ebx # load counter
        movl %1,%%ecx # load &x
        movl %2,%%edx # load &y
        clc # clear carry
%=:
        movl -4(%%edx,%%ebx,4),%%eax # load y[i]
        adcl %%eax,-4(%%ecx,%%ebx,4) # add with carry
        decl %%ebx # decrement counter
        jnz %=b # loop"
   :
   : "g"(N), "g"(x.data), "g"(y.data)
   : "eax", "ebx", "ecx", "edx", "memory", "cc");
  }

  template <std::size_t N>
  void kernel<N>::add(num_t& x, const uint32_t& y) {
    asm("movl %0,%%ebx # load counter
        movl %1,%%ecx # load x
        addl %2,(%%ecx,%%ebx,4) # add y
%=0:
        jnc %=1f # abort loop if carry not set
        adcl  $0,-4(%%ecx,%%ebx,4) # add with carry
        decl %%ebx # decrement counter
        jnz %=0b # loop
%=1:"
         :
         : "g"(N-1), "g"(x.data), "r"(y)
         : "ebx", "ecx", "memory", "cc");
  }

  template <std::size_t N>
  void kernel<N>::add(num_t& x, const uint64_t& y) {
    uint64_hilo tmp;
    tmp.i64 = y;
    asm("movl %0,%%ebx # load counter
        movl %1,%%ecx # load x
        addl %2,4(%%ecx,%%ebx,4)
        adcl %3,(%%ecx,%%ebx,4)
%=0:
        jnc %=1f # abort loop if carry not set
        adcl  $0,-4(%%ecx,%%ebx,4)
        decl %%ebx # decrement counter
        jnz %=0b # loop
%=1:"
         :
         : "g"(N-2), "g"(x.data), "r"(tmp.i32[0]), "r"(tmp.i32[1])
         : "ebx", "ecx", "memory", "cc");
  }

  template <std::size_t N>
  void kernel<N>::sub(num_t& x, const num_t& y) {
    asm("movl %0,%%ebx # load counter
        movl %1,%%ecx # load x
        movl %2,%%edx # load y
        clc # clear carry
%=:
        movl -4(%%edx,%%ebx,4),%%eax
        sbbl %%eax,-4(%%ecx,%%ebx,4)
        decl %%ebx # decrement counter
        jnz %=b # loop"
   :
   : "g"(N), "g"(x.data), "g"(y.data)
   : "eax", "ebx", "ecx", "edx", "memory", "cc");
  }

  template <std::size_t N>
  void kernel<N>::sub(num_t& x, const uint32_t& y) {
    asm("movl %0,%%ebx # load counter
        movl %1,%%ecx # load x
        subl %2,(%%ecx,%%ebx,4)
%=0:
        jnc %=1f # abort loop if carry not set
        sbbl  $0,-4(%%ecx,%%ebx,4)
        decl %%ebx # decrement counter
        jnz %=0b # loop
%=1:"
         :
         : "g"(N-1), "g"(x.data), "r"(y)
         : "ebx", "ecx", "memory", "cc");
  }

  template <std::size_t N>
  void kernel<N>::sub(num_t& x, const uint64_t& y) {
    uint64_hilo tmp;
    tmp.i64 = y;
    asm("movl %0,%%ebx # load counter
        movl %1,%%ecx # load x
        subl %2,4(%%ecx,%%ebx,4)
        sbbl %3,(%%ecx,%%ebx,4)
%=0:
        jnc %=1f # abort loop if carry not set
        sbbl  $0,-4(%%ecx,%%ebx,4)
        decl %%ebx # decrement counter
        jnz %=0b # loop
%=1:"
         :
         : "g"(N-2), "g"(x.data), "r"(tmp.i32[0]), "r"(tmp.i32[1])
         : "ebx", "ecx", "memory", "cc");
  }

  template <std::size_t N>
  void kernel<N>::mul(num_t& w, const num_t& x, const num_t& y) {
    init(w);
    for ( int j = N-1; j >= 0; --j ) {
      uint32_t carry = 0UL;
      for ( int i = N-1; i >= 0 && i+j-N+1 >= 0; --i ) {
  uint64_hilo tmp;
  tmp.i64 = static_cast<uint64_t>(x.data[j]) * y.data[i]
    + static_cast<uint64_t>(w.data[i+j-N+1]) + carry;
        w.data[i+j-N+1] = tmp.i32[0];
        carry = tmp.i32[1];
      }
    }
  }

  template <std::size_t N>
  void kernel<N>::mul(num_t& w, const num_t& x, const uint32_t& y) {
    uint32_t carry = 0UL;
    for ( int i = N-1; i >= 0; --i ) {
      uint64_hilo tmp;
      tmp.i64 = static_cast<uint64_t>(y) * x.data[i] + carry;
      w.data[i] = tmp.i32[0];
      carry = tmp.i32[1];
    }
  }

  template <std::size_t N>
  void kernel<N>::mul(num_t& w, const num_t& x, const uint64_t& y) {
    uint64_hilo tmp1;
    tmp1.i64 = y;
    const uint32_t ylo = tmp1.i32[0];
    const uint32_t yhi = tmp1.i32[1];
    uint32_t carry = 0UL;
    for ( int i = N-1; i >= 0; --i ) {
      uint64_hilo tmp;
      tmp.i64 = static_cast<uint64_t>(ylo) * x.data[i] + carry;
      w.data[i] = tmp.i32[0];
      carry = tmp.i32[1];
    }
    carry = 0UL;
    for ( int i = N-1; i > 0; --i ) {
      uint64_hilo tmp;
      tmp.i64 = static_cast<uint64_t>(yhi) * x.data[i]
  + static_cast<uint64_t>(w.data[i-1]) + carry;
      w.data[i-1] = tmp.i32[0];
      carry = tmp.i32[1];
    }
  }

  template <std::size_t N>
  void kernel<N>::div(num_t& q, uint32_t& r, const num_t& x, const uint32_t& y) {
    r = 0UL;
    for ( int i = 0; i < N ; ++i ) {
      uint64_t tmp = (static_cast<uint64_t>(r) << 32) + x.data[i];
      q.data[i] = tmp / y;
      r = tmp % y;
    }
  }

  template <std::size_t N>
  bool kernel<N>::less(const num_t& x, const num_t& y) {
    for ( int i = 0; i < N; ++i ) {
      if ( x.data[i] > y.data[i] )
  return false;
      if ( x.data[i] < y.data[i] )
  return true;
    }
    return false;
  }

  template <std::size_t N>
  bool kernel<N>::less(const num_t& x, const uint32_t& y) {
    for ( int i = 0; i < N-1; ++i )
      if ( x.data[i] != 0UL )
  return false;
    if ( x.data[N-1] < y )
      return true;
    return false;
  }

  template <std::size_t N>
  bool kernel<N>::less(const num_t& x, const uint64_t& y) {
    for ( int i = 0; i < N-2; ++i )
      if ( x.data[i] != 0UL )
  return false;
    uint64_hilo tmp;
    tmp.i64 = y;
    if ( x.data[N-2] < tmp.i32[1] )
      return true;
    if ( x.data[N-2] > tmp.i32[1] )
      return false;
    if ( x.data[N-1] < tmp.i32[0] )
      return true;
    return false;
  }

  template <std::size_t N>
  bool kernel<N>::equal(const num_t& x, const num_t& y) {
    for ( int i = 0; i < N; ++i )
      if ( x.data[i] != y.data[i] )
  return false;
    return true;
  }

  template <std::size_t N>
  bool kernel<N>::equal(const num_t& x, const uint32_t& y) {
    for ( int i = 0; i < N-1; ++i )
      if ( x.data[i] != 0UL )
  return false;
    if ( x.data[N-1] != y )
      return false;
    return true;
  }

  template <std::size_t N>
  bool kernel<N>::equal(const num_t& x, const uint64_t& y) {
    for ( int i = 0; i < N-2; ++i )
      if ( x.data[i] != 0UL )
  return false;
    uint64_hilo tmp;
    tmp.i64 = y;
    if ( x.data[N-2] != tmp.i32[1] )
      return false;
    if ( x.data[N-1] != tmp.i32[0] )
      return false;
    return true;
  }

  template <std::size_t N>
  bool kernel<N>::greater(const num_t& x, const num_t& y) {
    for ( int i = 0; i < N; ++i ) {
      if ( x.data[i] < y.data[i] )
  return false;
      if ( x.data[i] > y.data[i] )
  return true;
    }
    return false;
  }

  template <std::size_t N>
  bool kernel<N>::greater(const num_t& x, const uint32_t& y) {
    for ( int i = 0; i < N-1; ++i )
      if ( x.data[i] != 0UL )
  return true;
    if ( x.data[N-1] > y )
      return true;
    return false;
  }

  template <std::size_t N>
  bool kernel<N>::greater(const num_t& x, const uint64_t& y) {
    for ( int i = 0; i < N-2; ++i )
      if ( x.data[i] != 0UL )
  return true;
    uint64_hilo tmp;
    tmp.i64 = y;
    if ( x.data[N-2] > tmp.i32[1] )
      return true;
    if ( x.data[N-2] < tmp.i32[1] )
      return false;
    if ( x.data[N-1] > tmp.i32[0] )
      return true;
    return false;
  }

}

#endif
