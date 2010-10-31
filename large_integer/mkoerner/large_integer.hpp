#ifndef LARGE_INTEGER_HPP
#define LARGE_INTEGER_HPP

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

#include <algorithm>
#include <string>
#include <cstddef>

namespace symbolic {
  
  /*
   * The class lint<NLIMB,KERNEL> is a high precision
   * signed integer class. The internals (programmed in
   * assembler, C , C++ or whatever) are encapsulated in
   * the class  which (at some point) should have
   * a reasonable and fixed interface, so that it can be
   * replaced by different kernels on different machines.
   * the number NLIMB is the number of elementary storage
   * units (the KERNEL specifies, which type that is)
   * are allocated per lint.
   */

  template <std::size_t N, class KERNEL>
  class lint {
  public:
    
    typedef short int sign_type;
    typedef typename KERNEL::num_t num_t;
    
    // constructors
    lint();
    lint(const int&);
    lint(const long long&);
    template <std::size_t M, class KERNEL2>
    lint(const lint<M,KERNEL2>&);
    explicit lint(const std::string&);

    // operations with lints
    inline lint& operator+=(const lint&);
    inline lint& operator-=(const lint&);
    inline lint& operator*=(const lint&);

    inline lint& operator-() const;

    // operations with int
    inline lint& operator+=(const int&);
    inline lint& operator-=(const int&);
    inline lint& operator*=(const int&);

    inline lint& operator+=(const long long&);
    inline lint& operator-=(const long long&);
    inline lint& operator*=(const long long&);

    inline bool operator==(const lint&) const;
    inline bool operator!=(const lint&) const;

    inline bool operator<(const lint&) const;

    template <class Ostream>
    Ostream& print(Ostream&) const;

    //private:
    
    sign_type sign;
    num_t l;

  };

  template <std::size_t N, class L>
  std::ostream& operator<<(std::ostream&, const lint<N,L>&);

  template <std::size_t N, class L>
  std::istream& operator>>(std::istream&, lint<N,L>&);
  
  template <std::size_t N, class L>
  lint<N,L> operator+(const lint<N,L>&, const lint<N,L>&);
  template <std::size_t N, class L>
  lint<N,L> operator-(const lint<N,L>&, const lint<N,L>&);
  template <std::size_t N, class L>
  lint<N,L> operator*(const lint<N,L>&, const lint<N,L>&);
  
  template <std::size_t N, class L>
  lint<N,L> operator+(const lint<N,L>&, const int&);
  template <std::size_t N, class L>
  lint<N,L> operator-(const lint<N,L>&, const int&);
  template <std::size_t N, class L>
  lint<N,L> operator*(const lint<N,L>&, const int&);
  template <std::size_t N, class L>
  lint<N,L> operator+(const int&, const lint<N,L>&);
  template <std::size_t N, class L>
  lint<N,L> operator-(const int&, const lint<N,L>&);
  template <std::size_t N, class L>
  lint<N,L> operator*(const int&, const lint<N,L>&);

  template <std::size_t N, class L>
  lint<N,L> operator+(const lint<N,L>&, const long long&);
  template <std::size_t N, class L>
  lint<N,L> operator-(const lint<N,L>&, const long long&);
  template <std::size_t N, class L>
  lint<N,L> operator*(const lint<N,L>&, const long long&);
  template <std::size_t N, class L>
  lint<N,L> operator+(const long long&, const lint<N,L>&);
  template <std::size_t N, class L>
  lint<N,L> operator-(const long long&, const lint<N,L>&);
  template <std::size_t N, class L>
  lint<N,L> operator*(const long long&, const lint<N,L>&);

  /////////////////////////IMPLEMENTATION////////////////////////

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>::lint() {}

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>::lint(const int& x) {
    if ( x < 0 ) {
      sign = -1;
      KERNEL::init(l,-x);
    } else {
      sign = 1;
      KERNEL::init(l,x);
    }
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>::lint(const long long& x) {
    if ( x < 0 ) {
      sign = -1;
      KERNEL::init(l,-x);
    } else {
      sign = 1;
      KERNEL::init(l,x);
    }
  }

  template <std::size_t N, class KERNEL> template <std::size_t M, class KERNEL2>
  lint<N,KERNEL>::lint(const lint<M,KERNEL2>& x) {
    sign = x.sign;
    KERNEL::extend_init(l,x.l);
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>::lint(const std::string& t) {
    if ( t[0] == '-' ) {
      sign = -1;
      KERNEL::init(l,t.substr(1).c_str());
    } else {
      sign = 1;
      KERNEL::init(l,t.c_str());
    }
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator+=(const lint& x) {
    if ( sign == x.sign ) {
      KERNEL::add(l,x.l);
      return *this;
    }
    if ( KERNEL::less(l,x.l) ) {
      sign = x.sign;
      num_t t = l;
      l = x.l;
      KERNEL::sub(l,t);
      return *this;
    }
    KERNEL::sub(l,x.l);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator-=(const lint& x) {
    if ( sign != x.sign ) {
      KERNEL::add(l,x.l);
      return *this;
    }
    if ( KERNEL::less(l,x.l) ) {
      sign = -x.sign;
      num_t t = l;
      l = x.l;
      KERNEL::sub(l,t);
      return *this;
    }
    KERNEL::sub(l,x.l);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator*=(const lint& x) {
    sign *= x.sign;
    num_t t = l;
    KERNEL::mul(l,t,x.l);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator-() const {
    lint<N,KERNEL> tmp(*this);
    this->sign = -this->sign;
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator+=(const int& x) {
    int xsign, xabs;
    if ( x < 0 ) {
      xsign = -1;
      xabs = -x;
    } else {
      xsign = 1;
      xabs = x;
    }
    if ( sign == xsign ) {
      KERNEL::add(l,xabs);
      return *this;
    }
    if ( KERNEL::less(l,xabs) ) {
      sign = xsign;
      num_t tmp = l;
      KERNEL::init(l,xabs);
      KERNEL::sub(l,tmp);
      return *this;
    }
    KERNEL::sub(l,xabs);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator-=(const int& x) {
    int xsign, xabs;
    if ( x < 0 ) {
      xsign = -1;
      xabs = -x;
    } else {
      xsign = 1;
      xabs = x;
    }
    if ( sign != xsign ) {
      KERNEL::add(l,xabs);
      return *this;
    }
    if ( KERNEL::less(l,xabs) ) {
      sign = -xsign;
      num_t tmp = l;
      KERNEL::init(l,xabs);
      KERNEL::sub(l,tmp);
      return *this;
    }
    KERNEL::sub(l,xabs);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator*=(const int& x) {
    int xsign, xabs;
    if ( x < 0 ) {
      xsign = -1;
      xabs = -x;
    } else {
      xsign = 1;
      xabs = x;
    }
    sign *= xsign;
    num_t t = l;
    KERNEL::mul(l,t,xabs);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator+=(const long long& x) {
    int xsign;
    long long xabs;
    if ( x < 0 ) {
      xsign = -1;
      xabs = -x;
    } else {
      xsign = 1;
      xabs = x;
    }
    if ( sign == xsign ) {
      KERNEL::add(l,xabs);
      return *this;
    }
    if ( KERNEL::less(l,xabs) ) {
      sign = xsign;
      num_t tmp = l;
      KERNEL::init(l,xabs);
      KERNEL::sub(l,tmp);
      return *this;
    }
    KERNEL::sub(l,xabs);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator-=(const long long& x) {
    int xsign;
    long long xabs;
    if ( x < 0 ) {
      xsign = -1;
      xabs = -x;
    } else {
      xsign = 1;
      xabs = x;
    }
    if ( sign != xsign ) {
      KERNEL::add(l,xabs);
      return *this;
    }
    if ( KERNEL::less(l,xabs) ) {
      sign = -xsign;
      num_t tmp = l;
      KERNEL::init(l,xabs);
      KERNEL::sub(l,tmp);
      return *this;
    }
    KERNEL::sub(l,xabs);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  lint<N,KERNEL>& lint<N,KERNEL>::operator*=(const long long& x) {
    int xsign;
    long long xabs;
    if ( x < 0 ) {
      xsign = -1;
      xabs = -x;
    } else {
      xsign = 1;
      xabs = x;
    }
    sign *= xsign;
    num_t t = l;
    KERNEL::mul(l,t,xabs);
    return *this;
  }

  template <std::size_t N, class KERNEL>
  bool lint<N,KERNEL>::operator==(const lint& x) const {
    return KERNEL::equal(l,x.l);
  }

  template <std::size_t N, class KERNEL>
  bool lint<N,KERNEL>::operator!=(const lint& x) const {
    return ! ( *this == x );
  }

  template <std::size_t N, class KERNEL>
  bool lint<N,KERNEL>::operator<(const lint& x) const {
    // this routine is not entirely correct!
    if ( sign < x.sign )
      return true;
    if ( sign > x.sign )
      return false;
    if ( sign == 1 && KERNEL::less(l,x.l) )
      return true;
    if ( sign == -1 && KERNEL::less(l,x.l) )
      return true;
    return false;
  } 
  
  template <std::size_t N, class KERNEL> template <class Ostream>
  Ostream& lint<N,KERNEL>::print(Ostream& os) const {
    if ( sign == -1 )
      os << '-';
    KERNEL::print(l,os);
    return os;
  }
  
  template <std::size_t N, class L>
  std::ostream& operator<<(std::ostream& os, const lint<N,L>& l) {
    return l.print(os);
  }

  template <std::size_t N, class L>
  std::istream& operator>>(std::istream& is, lint<N,L>& l) {
    std::string tmp;
    is >> tmp;
    l = lint<N,L>(tmp);
    return is;
  }

  template <std::size_t N, class L>
  lint<N,L> operator+(const lint<N,L>& x, const lint<N,L>& y) {
    lint<N,L> temp(x);
    temp += y;
    return temp;
  }

  template <std::size_t N, class L>
  lint<N,L> operator-(const lint<N,L>& x, const lint<N,L>& y) {
    lint<N,L> temp(x);
    temp -= y;
    return temp;
  }
  
  template <std::size_t N, class L>
  lint<N,L> operator*(const lint<N,L>& x, const lint<N,L>& y) {
    lint<N,L> temp(x);
    temp *= y;
    return temp;
  }
  
  template <std::size_t N, class L>
  lint<N,L> operator+(const lint<N,L>& x, const int& y) {
    lint<N,L> temp(x);
    temp += y;
    return temp;
  }
  
  template <std::size_t N, class L>
  lint<N,L> operator-(const lint<N,L>& x, const int& y) {
    lint<N,L> temp(x);
    temp -= y;
    return temp;
  }
  
  template <std::size_t N, class L>
  lint<N,L> operator*(const lint<N,L>& x, const int& y) {
    lint<N,L> temp(x);
    temp *= y;
    return temp;
  }
  
  template <std::size_t N, class L>
  lint<N,L> operator+(const lint<N,L>& x, const long long& y) {
    lint<N,L> temp(x);
    temp += y;
    return temp;
  }
  
  template <std::size_t N, class L>
  lint<N,L> operator-(const lint<N,L>& x, const long long& y) {
    lint<N,L> temp(x);
    temp -= y;
    return temp;
  }
  
  template <std::size_t N, class L>
  lint<N,L> operator*(const lint<N,L>& x, const long long& y) {
    lint<N,L> temp(x);
    temp *= y;
    return temp;
  }

  template <std::size_t N, class L>
  lint<N,L> operator+(const int&, const lint<N,L>&);
  template <std::size_t N, class L>
  lint<N,L> operator-(const int&, const lint<N,L>&);
  template <std::size_t N, class L>
  lint<N,L> operator*(const int&, const lint<N,L>&);

}

#endif
