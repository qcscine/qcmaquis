#include <iostream>
#include "polynomial.hpp"

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

int main() {
  typedef symbolic::polynomial<int,'t',10> pf;
  typedef pf::polynom_type polynom;
  typedef pf::symbol_type symbol;
  typedef pf::value_type coefficient;

  symbol t;

  using std::cout;
  using std::endl;

  cout << t*t << endl;
  cout << t+t << endl;
  cout << t-t << endl;
  cout << 1+t << endl;
  cout << t+1 << endl;
  cout << 1-t << endl;
  cout << t-1 << endl;
  cout << 1*t << endl;
  cout << t*1 << endl;
  cout << (t==t) << endl;
  cout << (t!=t) << endl;
  cout << ((1+t)*(1-t)*t*t==t*t) << endl;
  cout << ((1+t)*(1-t)*t*t!=t*t) << endl;

  polynom p(1+t*t);

  cout << t+p << endl;
  cout << t-p << endl;
  cout << t*p << endl;
  cout << (p+=1) << endl;
  cout << p+1 << endl;
  cout << 1+p << endl;
  cout << (p-=1) << endl;
  cout << p-1 << endl;
  cout << 1-p << endl;
  cout << (p*=1) << endl;
  cout << p*1 << endl;
  cout << 1*p << endl;
  cout << p*t << endl;
  cout << t*p << endl;
  cout << 1+p << endl;
  cout << 1-p << endl;

  return 0;
}
