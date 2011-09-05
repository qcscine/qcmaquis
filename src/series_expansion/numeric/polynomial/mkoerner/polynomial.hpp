#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

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

/****************************************************************
 *
 * This is a simple class for polynomial arithmetic in one
 * variable. It can be used as follows:
 *
 * example.C:
 *
 * typedef symbolic::polynomial<int,'t',10> pf;
 * typedef pf::polynom_type polynom;
 * typedef pf::symbol_type symbol;
 * typedef pf::value_type coefficient;           
 * symbol t;
 * polynom p(t*t);
 * ...
 *
 * We can separate the occuring types as follows:
 * coefficients, monoms and polynoms.
 * coefficients are simple numbers. The type of the
 * coefficients is determined by the template parameter
 * CoefficientType. A monom is an expression of the form
 * c*x^n where x is a symbol. A symbol is simply represented
 * as a monom of degeree 1. The monoms have to know the maximal
 * degree of the resulting polynom for addition and subtraction.
 * This is also the reason for the creation of types using
 * polynomial<CoefficientType, Symbol, MaxDegree>.
 * The polynom class represents does not check for bounds
 * violations. Note however that an expression like
 * polynom<...,md>+monom<...,md,n> with n>md will result in
 * a segmentation fault.
 *
 * (c) 2000 Mathias Koerner
 *
 ****************************************************************/

#include<iosfwd>

namespace symbolic {

  template<class ct, char s, int md, int n> class monom;
  template<class ct, char s, int md> class polynom;

  template<class CoefficientType,
    char Symbol,
    int MaxDegree>
  struct polynomial{
    typedef polynom<CoefficientType,
      Symbol,
      MaxDegree> polynom_type;
    typedef monom<CoefficientType,
      Symbol,
      MaxDegree,1> symbol_type;
    typedef CoefficientType value_type;
  };

}

template<class ct, char s, int md, int n>
std::ostream& operator<<(std::ostream&, const symbolic::monom<ct,s,md,n>&);

template<class ct, char s, int md>
std::ostream& operator<<(std::ostream&, const symbolic::polynom<ct,s,md>&);

namespace symbolic {

  template<class ct, char s, int md, int n> 
  class monom {
  public:
    monom()
      : c(static_cast<ct>(1)) {};
    explicit monom(const ct&);

    bool operator==(const polynom<ct,s,md>&) const;
    bool operator!=(const polynom<ct,s,md>&) const;

    template<int m>
    polynom<ct,s,md> operator+(const monom<ct,s,md,m>&) const;
    polynom<ct,s,md> operator+(const polynom<ct,s,md>&) const;
    polynom<ct,s,md> operator+(const ct&) const;

    template<int m>
    polynom<ct,s,md> operator-(const monom<ct,s,md,m>&) const;
    polynom<ct,s,md> operator-(const polynom<ct,s,md>&) const;
    polynom<ct,s,md> operator-(const ct&) const;

    template<int m>
    monom<ct,s,md,n+m> operator*(const monom<ct,s,md,m>&) const;
    polynom<ct,s,md> operator*(const polynom<ct,s,md>&) const;
    monom<ct,s,md,n> operator*(const ct&) const;

    ct& coefficient();
    const ct& coefficient() const;

  private:
    ct c;
  };

  template<class ct, char s, int md, int n>
  polynom<ct,s,md> operator+(const ct&, const monom<ct,s,md,n>&);

  template<class ct, char s, int md, int n>
  polynom<ct,s,md> operator-(const ct&, const monom<ct,s,md,n>&);

  template<class ct, char s, int md, int n>
  monom<ct,s,md,n> operator*(const ct&, const monom<ct,s,md,n>&);

  template<class ct, char s, int md>
  class polynom {
  public:
    polynom();
    explicit polynom(const ct&);
    template<int n>
    polynom(const monom<ct,s,md,n>&);

    polynom<ct,s,md>& operator=(const polynom<ct,s,md>&);

    bool operator==(const ct&) const;
    template<int n>
    bool operator==(const monom<ct,s,md,n>&) const;
    bool operator==(const polynom<ct,s,md>&) const;

    bool operator!=(const ct&) const;
    template<int n>
    bool operator!=(const monom<ct,s,md,n>&) const;
    bool operator!=(const polynom<ct,s,md>&) const;

    template<int n>
    polynom<ct,s,md>& operator+=(const monom<ct,s,md,n>&);
    polynom<ct,s,md>& operator+=(const polynom<ct,s,md>&);
    polynom<ct,s,md>& operator+=(const ct&);
    template<int n>
    polynom<ct,s,md> operator+(const monom<ct,s,md,n>&) const;
    polynom<ct,s,md> operator+(const polynom<ct,s,md>&) const;
    polynom<ct,s,md> operator+(const ct&) const;

    template<int n>
    polynom<ct,s,md>& operator-=(const monom<ct,s,md,n>&);
    polynom<ct,s,md>& operator-=(const polynom<ct,s,md>&);
    polynom<ct,s,md>& operator-=(const ct&);
    template<int n>
    polynom<ct,s,md> operator-(const monom<ct,s,md,n>&) const;
    polynom<ct,s,md> operator-(const polynom<ct,s,md>&) const;
    polynom<ct,s,md> operator-(const ct&) const;
   
    template<int n>
    polynom<ct,s,md> operator*(const monom<ct,s,md,n>&) const;
    polynom<ct,s,md> operator*(const polynom<ct,s,md>&) const;
    polynom<ct,s,md> operator*(const ct&) const;
    template<int n>
    polynom<ct,s,md>& operator*=(const monom<ct,s,md,n>&);
    polynom<ct,s,md>& operator*=(const polynom<ct,s,md>&);
    polynom<ct,s,md>& operator*=(const ct&);

    polynom<ct,s,md> operator-() const;

    ct& operator[](int i);
    const ct& operator[](int i) const;

    friend std::ostream& ::operator<<(std::ostream&, const polynom<ct,s,md>&); 
  private:
    ct c[md+1];
  };

  template<class ct, char s, int md>
  polynom<ct,s,md> operator+(const ct&, const polynom<ct,s,md>&);

  template<class ct, char s, int md>
  polynom<ct,s,md> operator-(const ct&, const polynom<ct,s,md>&);

  template<class ct, char s, int md>
  polynom<ct,s,md> operator*(const ct&, const polynom<ct,s,md>&);

}

#include "polynomial.cpp"

#endif
