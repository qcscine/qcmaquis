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

namespace symbolic {

  template<class ct, char s, int md, int n> 
  monom<ct,s,md,n>::monom(const ct& v)
    : c(v) {}

  template<class ct, char s, int md, int n>
  bool monom<ct,s,md,n>::operator==(const polynom<ct,s,md>& rhs) const {
    return static_cast<polynom<ct,s,md> >(*this)==rhs;
  }

  template<class ct, char s, int md, int n>
  bool monom<ct,s,md,n>::operator!=(const polynom<ct,s,md>& rhs) const {
    return !((*this)==rhs);
  }

  template<class ct, char s, int md, int n> template<int m>
  polynom<ct,s,md> monom<ct,s,md,n>::operator+(const monom<ct,s,md,m>& rhs) const {
    polynom<ct,s,md> buf;
    buf[n]+=c;
    buf[m]+=rhs.coefficient();
    return buf;
  }

  template<class ct, char s, int md, int n>
  polynom<ct,s,md> monom<ct,s,md,n>::operator+(const polynom<ct,s,md>& rhs) const {
    polynom<ct,s,md> buf(rhs);
    buf[n]+=c;
    return buf;
  }

  template<class ct, char s, int md, int n>
  polynom<ct,s,md> monom<ct,s,md,n>::operator+(const ct& rhs) const {
    return (*this)+static_cast<monom<ct,s,md,0> >(rhs); 
  }

  template<class ct, char s, int md, int n>
  polynom<ct,s,md> operator+(const ct& lhs, const monom<ct,s,md,n>& rhs) {
    return static_cast<monom<ct,s,md,0> >(lhs)+rhs;
  }

  template<class ct, char s, int md, int n> template<int m>
  polynom<ct,s,md> monom<ct,s,md,n>::operator-(const monom<ct,s,md,m>& rhs) const {
    polynom<ct,s,md> buf;
    buf[n]+=c;
    buf[m]-=rhs.coefficient();
    return buf;
  }

  template<class ct, char s, int md, int n>
  monom<ct,s,md,n> monom<ct,s,md,n>::operator-() const {
    return static_cast<monom<ct,s,md,n> >(0)-(*this);
  }

  template<class ct, char s, int md, int n>
  polynom<ct,s,md> monom<ct,s,md,n>::operator-(const polynom<ct,s,md>& rhs) const {
    polynom<ct,s,md> buf(-rhs);
    buf[n]+=c;
    return buf;
  }

  template<class ct, char s, int md, int n>
  polynom<ct,s,md> monom<ct,s,md,n>::operator-(const ct& rhs) const {
    return (*this)-static_cast<monom<ct,s,md,0> >(rhs); 
  }

  template<class ct, char s, int md, int n>
  polynom<ct,s,md> operator-(const ct& lhs, const monom<ct,s,md,n>& rhs) {
    return static_cast<monom<ct,s,md,0> >(lhs)-rhs;
  }

  template<class ct, char s, int md, int n> template<int m>
  monom<ct,s,md,n+m> monom<ct,s,md,n>::operator*(const monom<ct,s,md,m>& rhs) const {
    return monom<ct,s,md,n+m>(c*rhs.coefficient());
  }

  template<class ct, char s, int md, int n> 
  polynom<ct,s,md> monom<ct,s,md,n>::operator*(const polynom<ct,s,md>& rhs) const {
    polynom<ct,s,md> buf;
    if (!rhs.is_zero()) {
      for(int i=0; i<=md-n; ++i)
	buf[i+n]=c*rhs[i];
      for(int i=0; i<n; ++i)
	buf[i]=static_cast<ct>(0);
    }
    return buf;
  }

  template<class ct, char s, int md, int n>
  monom<ct,s,md,n> monom<ct,s,md,n>::operator*(const ct& rhs) const {
    return (*this)*static_cast<monom<ct,s,md,0> >(rhs); 
  }

  template<class ct, char s, int md, int n>
  monom<ct,s,md,n> operator*(const ct& lhs, const monom<ct,s,md,n>& rhs) {
    return static_cast<monom<ct,s,md,0> >(lhs)*rhs;
  }

  template<class ct, char s, int md, int n>
  ct& symbolic::monom<ct,s,md,n>::coefficient() {
    return c;
  }
  
  template<class ct, char s, int md, int n>
  const ct& symbolic::monom<ct,s,md,n>::coefficient() const {
    return c;
  }

  template<class ct, char s, int md>
  polynom<ct,s,md>::polynom()
    : zero(true) {}

  template<class ct, char s, int md>
  polynom<ct,s,md>::polynom(const ct& v) {
    if (v==0)
      zero=true;
    else {
      c[0]=v;
      for (int i=1; i<=md; ++i)
	c[i]=static_cast<ct>(0);
      zero=false;
    }
  }

  template<class ct, char s, int md> template <int n>
  polynom<ct,s,md>::polynom(const monom<ct,s,md,n>& rhs) {
    if (rhs.coefficient()==static_cast<ct>(0))
      zero=true;
    else {
      invalidate_zero();
      c[n]=rhs.coefficient();
    }
  }

  template<class ct, char s, int md>
  polynom<ct,s,md>& polynom<ct,s,md>::operator=(const polynom<ct,s,md>& rhs) {
    if (rhs.is_zero())
      zero=true;
    else {
      for (int i=0; i<=md; ++i)
        c[i]=rhs.c[i];
      zero=false;
    }
    return *this;
  }

  template<class ct, char s, int md>
  bool polynom<ct,s,md>::operator==(const ct& rhs) const {
    if (zero && (rhs==0))
      return true;
    if (c[0]!=rhs)
      return false;
    for (int i=1; i<md; ++i)
      if (c[i]!=static_cast<ct>(0))
	return false;
    return true;
  }

  template<class ct, char s, int md> template<int n>
  bool polynom<ct,s,md>::operator==(const monom<ct,s,md,n>& rhs) const {
    if(c[n]!=rhs.coefficient())
      return false;
    for (int i=0; i<n; ++i)
      if (c[i]!=static_cast<ct>(0))
	return false;
    for (int i=n+1; i<md; ++i)
      if (c[i]!=static_cast<ct>(0))
	return false;
    return true;
  }

  template<class ct, char s, int md> 
  bool polynom<ct,s,md>::operator==(const polynom<ct,s,md>& rhs) const {
    if (zero && rhs.zero)
      return true;
    for (int i=0; i<=md; ++i)
      if (c[i]!=rhs[i])
	return false;
    return true;
  }

  template<class ct, char s, int md> 
  bool polynom<ct,s,md>::operator!=(const ct& rhs) const {
    return !((*this)==rhs);
  }

  template<class ct, char s, int md> template<int n>
  bool polynom<ct,s,md>::operator!=(const monom<ct,s,md,n>& rhs) const {
    return !((*this)==rhs);
  }

  template<class ct, char s, int md> 
  bool polynom<ct,s,md>::operator!=(const polynom<ct,s,md>& rhs) const {
    return !((*this)==rhs);
  }

  template<class ct, char s, int md> template<int n>
  polynom<ct,s,md>& polynom<ct,s,md>::operator+=(const monom<ct,s,md,n>& rhs) {
    if (rhs.coefficient()!=static_cast<ct>(0)) {
      if (zero)
	zero_invalidate();
      c[n]+=rhs.coefficient();
    }
    return *this;
  }

  template<class ct, char s, int md> 
  polynom<ct,s,md>& polynom<ct,s,md>::operator+=(const polynom<ct,s,md>& rhs) {
    if (!rhs.zero) {
      if (zero)
	zero_invalidate();
      for (int i=0; i<=md; ++i)
	c[i]+=rhs.c[i];
    }
    return *this;
  }

  template<class ct, char s, int md> 
  polynom<ct,s,md>& polynom<ct,s,md>::operator+=(const ct& rhs) {
    return (*this)+=static_cast<monom<ct,s,md,0> >(rhs);
  }
  
  template<class ct, char s, int md> template<int n>
  polynom<ct,s,md> polynom<ct,s,md>::operator+(const monom<ct,s,md,n>& rhs) const {
    polynom<ct,s,md> buf(*this);
    buf+=rhs;
    return buf;
  }

  template<class ct, char s, int md> 
  polynom<ct,s,md> polynom<ct,s,md>::operator+(const polynom<ct,s,md>& rhs) const {
    polynom<ct,s,md> buf(*this);
    buf+=rhs;
    return buf;
  }

  template<class ct, char s, int md>
  polynom<ct,s,md> polynom<ct,s,md>::operator+(const ct& rhs) const {
    polynom<ct,s,md> buf(*this);
    buf+=rhs;
    return buf;
  }

  template<class ct, char s, int md>
  polynom<ct,s,md> operator+(const ct& lhs, const polynom<ct,s,md>& rhs) {
    return rhs+lhs;
  }

  template<class ct, char s, int md> template<int n>
  polynom<ct,s,md>& polynom<ct,s,md>::operator-=(const monom<ct,s,md,n>& rhs) {
    (*this)+=(-rhs);
    return *this;
  }

  template<class ct, char s, int md> 
  polynom<ct,s,md>& polynom<ct,s,md>::operator-=(const polynom<ct,s,md>& rhs) {
    if (!rhs.zero) {
      if (zero)
	zero_invalidate();
      for (int i=0; i<=md; ++i)
	c[i]-=rhs.c[i];
    }
    return *this;
  }
  
  template<class ct, char s, int md> 
  polynom<ct,s,md>& polynom<ct,s,md>::operator-=(const ct& rhs) {
    return (*this)-=static_cast<monom<ct,s,md,0> >(rhs);
  }
  
  template<class ct, char s, int md> template<int n>
  polynom<ct,s,md> polynom<ct,s,md>::operator-(const monom<ct,s,md,n>& rhs) const {
    polynom<ct,s,md> buf(*this);
    buf-=rhs;
    return buf;
  }

  template<class ct, char s, int md>
  polynom<ct,s,md> polynom<ct,s,md>::operator-(const polynom<ct,s,md>& rhs) const {
    polynom<ct,s,md> buf(*this);
    buf-=rhs;
    return buf;
  }

  template<class ct, char s, int md>
  polynom<ct,s,md> polynom<ct,s,md>::operator-(const ct& rhs) const {
    polynom<ct,s,md> buf(*this);
    buf-=rhs;
    return buf;
  }

  template<class ct, char s, int md>
  polynom<ct,s,md> operator-(const ct& lhs, const polynom<ct,s,md>& rhs) {
    polynom<ct,s,md> buf(-rhs);
    buf+=lhs;
    return buf;
  }

  template<class ct, char s, int md> template<int n>
  polynom<ct,s,md> polynom<ct,s,md>::operator*(const monom<ct,s,md,n>& rhs) const {
    polynom<ct,s,md> buf;
    if (rhs.coefficient()==static_cast<ct>(0))
      zero=true;
    else
      if (!zero) {
	for(int i=0; i<=md-n; ++i)
	  buf[i+n]=rhs.coefficient()*c[i];
	for(int i=0; i<n; ++i)
	  buf[i]=static_cast<ct>(0);
      }
    return buf;
  }

  template<class ct, char s, int md> 
  polynom<ct,s,md> polynom<ct,s,md>::operator*(const polynom<ct,s,md>& rhs) const {
    polynom<ct,s,md> buf;
    if ( (!zero) && (!rhs.zero) )
      for (int i=0; i<=md; ++i) {
	buf[i]=static_cast<ct>(0);
	for (int j=0; j<=i; ++j)
	  buf[i]+=c[j]*rhs.c[i-j];
      }
    return buf;
  }

  template<class ct, char s, int md> 
  polynom<ct,s,md> polynom<ct,s,md>::operator*(const ct& rhs) const {
    polynom<ct,s,md> buf(*this);
    if (!zero)
      for (int i=0; i<=md; ++i)
	buf[i]*=rhs;
    return buf;
  }

  template<class ct, char s, int md> template<int n>
  polynom<ct,s,md>& polynom<ct,s,md>::operator*=(const monom<ct,s,md,n>& rhs) {
    polynom<ct,s,md> buf((*this)*rhs);
    (*this)=buf;
    return *this;
  }

  template<class ct, char s, int md>
  polynom<ct,s,md>& polynom<ct,s,md>::operator*=(const polynom<ct,s,md>& rhs) {
    polynom<ct,s,md> buf((*this)*rhs);
    (*this)=buf;
    return *this;
  }

  template<class ct, char s, int md>
  polynom<ct,s,md>& polynom<ct,s,md>::operator*=(const ct& rhs) {
    polynom<ct,s,md> buf((*this)*rhs);
    (*this)=buf;
    return *this;
  }

  template<class ct, char s, int md> 
  polynom<ct,s,md> operator*(const ct& lhs, const polynom<ct,s,md>& rhs) {
    return rhs*lhs;
  }
  
  template<class ct, char s, int md> 
  polynom<ct,s,md> polynom<ct,s,md>::operator-() const {
    polynom<ct,s,md> buf;
    if (!zero)
      for (int i=0; i<=md; ++i)
	buf.c[i]=-c[i];
    return buf;
  }

  template<class ct, char s, int md>
  typename polynom<ct,s,md>::sprod_type polynom<ct,s,md>::ext_mult(const polynom<ct,s,md>& rhs) const {
    sprod_type buf;
    if ( (!zero) && (!rhs.zero) )
      for (int i=0; i<=2*md; ++i) {
	buf[i]=static_cast<ct>(0);
	for (int j=0; j<=i; ++j)
	  buf[i]+=c[j]*rhs.c[i-j];
      }
    return buf;
  }

  template<class ct, char s, int md>
  void polynom<ct,s,md>::set_zero() {
    zero=true;
  }

  template<class ct, char s, int md>
  bool polynom<ct,s,md>::is_zero() const {
    return zero;
  }

  template<class ct, char s, int md> 
  ct& polynom<ct,s,md>::operator[](int i) {
    zero_invalidate();
    return c[i];
  }

  template<class ct, char s, int md> 
  const ct polynom<ct,s,md>::operator[](int i) const {
    if (zero)
      return 0;
    return c[i];
  }

  template<class ct, char s, int md> 
  void polynom<ct,s,md>::zero_invalidate() {
    zero=false;
    set_array_zero();
  }

  template<class ct, char s, int md>
  void polynom<ct,s,md>::set_array_zero() {
    for (int i=0; i<=md; ++i)
      c[i]=0;
  }

}

template<class ct, char s, int md, int n>
std::ostream& operator<<(std::ostream& os, 
			 const symbolic::monom<ct,s,md,n>& m) {
  os << m.coefficient() << "*" << s << "^" << n;
  return os;
}

template<class ct, char s, int md>
std::ostream& operator<<(std::ostream& os, 
			 const symbolic::monom<ct,s,md,1>& m) {
  os << m.coefficient() << "*" << s;
  return os;
}

template<class ct, char s, int md>
std::ostream& operator<<(std::ostream& os, 
			 const symbolic::monom<ct,s,md,0>& m) {
  os << m.coefficient();
  return os;
}

template<class ct, char s, int md>
std::ostream& operator<<(std::ostream& os, 
			 const symbolic::polynom<ct,s,md>& p) {
  bool printedStuff=false;
  bool noStar=false;
  for (int i=md; i>=0; --i)
    if (p[i] != static_cast<ct>(0)) {
      if (!printedStuff) {
	if (p[i]==static_cast<ct>(1) && i>0)
	  noStar=true;
	else
	  if (p[i]==static_cast<ct>(-1) && i>0) {
	    noStar=true;
	    os << "-";
	  } else
	    os << p[i];
	printedStuff=true;
      } else {
	if (p[i] < static_cast<ct>(0))
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
