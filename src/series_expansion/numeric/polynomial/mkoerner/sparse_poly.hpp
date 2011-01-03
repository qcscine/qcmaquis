#ifndef POLYNOMIAL_OPT_HPP
#define POLYNOMIAL_OPT_HPP

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

#include <iostream>

namespace symbolic {
  
  template <typename T, char S, int N, int O> struct monom;
  template <typename T, char S, int N> class polynom;
  
  template<typename CoefficientType,
    char Symbol,
    int MaxEntries>
  struct polynomial {
    typedef polynom<CoefficientType,
      Symbol,
      MaxEntries> polynom_type;
    typedef monom<CoefficientType,
      Symbol,
      MaxEntries,1> symbol_type;
    typedef CoefficientType value_type;
  };

  template <typename T, char S, int N, int O> 
  struct monom {
  public:
    explicit monom(const T& x = 1) 
      : value(x) {}
    inline monom<T,S,N,O> operator*(const T& x) const {
      monom<T,S,N,O> buf;
      buf.value = this->value * x;
      return buf;
    }
    T value;
  };
  
  template <typename T, char S, int N>
  class polynom {
  public:

    polynom()
      : zero(true) {}

    explicit polynom(const T& x) {
      if ( x == 0 ) {
	zero = true;
      } else {
	zero = false;
	offset = 0;
	c[0] = v;
	for ( int i = 1; i < N; ++i )
	  c[i] = 0;
      }
    }	

    template <typename U, int M>
    explicit polynom(const polynom<U,S,M>& x) {
      assert( M <= N );
      zero = x.zero;
      if ( !zero ) {
	offset = x.offset;
	for ( int i = 0; i < M; ++i )
	  c[i] = x.c[i];
	for ( int i = M; i < N; ++i )
	  c[i] = 0;
      }
    }

    template <typename U, int M>
    inline polynom<T,S,N>& operator+=(const polynom<U,S,M>& x) {
      assert( M <= N );
      if ( !rhs.zero ) {
	if ( zero ) {
	  offset = x.offset;
	  for ( int i = 0; i < M; ++i )
	    c[i] = x.c[i];
	  for ( int i = M; i < N; ++i )
	    c[i] = 0;
	} else {
	  assert( offset == x.offset );
	  for ( int i = 0; i < M; ++i )
	    c[i] += x.c[i];
	}
      }
      return *this;
    }

    template <typename U, int M, typename V, int L> inline void 
    add_product(const polynom<U,S,M>& x, const monom<V,S,L,1>& m) {
      assert( M <= N );
      if ( !x.zero )
	if ( zero ) {
	  zero = false;
	  if ( x.offset == 0 ) {
	    offset = 1;
	    for ( int i = 0; i < M; ++i )
	      c[i] = x.c[i] * m.value;
	    for ( int i = M; i < N; ++i )
	      c[i] = 0;
	  } else {
	    assert( M < N );
	    offset = 0;
	    c[0] = 0;
	    for ( int i = 0; i < M; ++i )
	      c[i+1] = x.c[i] * m.value;
	    for ( int i = M; i < N-1; ++i )
	      c[i+1] = 0;
	  } // end of section for this == zero
	} else {
	  if ( x.offset == 0 ) {
	    assert( offset == 1);
	    for ( int i = 0; i < M; ++i )
	      c[i] += x.c[i] * m.value;
	  } else {
	    assert( offset == 0 );
	    assert( M < N );
	    for ( int i = 0; i < M; ++i )
	      c[i+1] += x.c[i] * m.value;
	  }
        } // end of section for this != zero
    }
  
    template <typename U, int M, typename V> inline void 
    add_product(const polynom<T,S,M>& x, const V& v) {
      assert( M <= N );
      if ( !x.zero )
	if ( zero ) {
	  zero = false;
	  offset = x.offset;
	  for ( int i = 0; i < M; ++i )
	    c[i] = x.c[i] * v;
	  for ( int i = M; i < N; ++i )
	    c[i] = 0;
	} else {
	  assert( offset == x.offset );
	  for ( int i = 0; i < M; ++i )
	    c[i] += x.c[i] * v;
	}
    }
  
    template <typename U, int M, typename V, int L> inline void 
    ext_mult_add(const polynom<U,S,M>& x, const polynom<V,S,L>& y) {
      assert( (M+L-1)*2 <= N*2-1 );
      if ( (!x.zero) && !(y.zero) ) {
	if ( zero ) {
	  for ( int i = 0; i < N; ++i )
	    c[i] = 0;
	}
	int localoff;
	if ( x.offset + y.offset == 2 ) {
	  offset = 0;
	  localoff = 1;
	} else {
	  offset =  x.offset + y.offset;
	  localoff = 0;
	}
	for ( int i = 0; i < M; ++i )
	  if ( x.c[i] != 0 )
	    for ( int j = 0; j < L; ++j )
	      if ( y.c[j] != 0 ) {
		assert( i+j+localoff < N ); //debugging only!!
		c[i+j+localoff] += x.c[i] * y.c[j];
	      }
      }
    }

    inline bool is_zero() const 
    { return zero; }
    
    // reference or value? value easier.
    inline T operator[](int i) const {
      if ( zero || (i%2) != offset ) {
	return 0;
      } else {
	assert( i/2 < N );
	return c[i/2];
      }
    }
    
    template <class OSTREAM>
    OSTREAM& print(OSTREAM&) const;
    
  private:
    bool zero;
    int offset;
    ct c[md];
    
  };

  template <typename T, char S, int N> template <class OSTREAM>
  OSTREAM& polynom<T,S,N>::print(OSTREAM& os) const {
    bool printedStuff=false;
    bool noStar=false;
    for ( int i = N*2-1; i >= 0; --i )
      if ( p[i] != static_cast<T>(0) ) {
	if ( !printedStuff ) {
	  if ( p[i] == static_cast<ct>(1) && i>0 )
	    noStar = true;
	  else
	    if ( p[i] == static_cast<ct>(-1) && i>0 ) {
	      noStar = true;
	      os << "-";
	    } else
	      os << p[i];
	  printedStuff = true;
	} else {
	  if ( p[i] < static_cast<ct>(0) )
	    if ( (p[i] == static_cast<ct>(-1)) &&
		 i>0 ) {
	      os << "-";
	      noStar=true;
	    } else
	      os << p[i];
	  else
	    if ( (p[i] == static_cast<ct>(1)) &&
		 i>0 ) {
	      os << "+";
	      noStar=true;
	    } else
	      os << "+" << p[i];
	}
	switch (i) {
	case 0:
	  break;
	case 1:
	  if (!noStar)
	    os << "*";
	  else
	    noStar=false;
	  os << s;
	  break;
	default:
	  if (!noStar)
	    os << "*";
	  else
	    noStar=false;
	  os << s << "^" << i;
	}
      }
    if (printedStuff==false)
      os << static_cast<ct>(0);
    return os;
  }

  template <typename T, char S, int N>
  std::ostream& operator(std::ostream& os, const polynom<T,S,N>& p) {
    return p.print(os);
  }

}

#endif
