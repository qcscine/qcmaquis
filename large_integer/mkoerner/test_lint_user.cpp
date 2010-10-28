#include <iostream>
#include <string>
#include "large_integer.hpp"
//#include "kernel_generic.hpp"
#include "kernel_i386.hpp"

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

const int N = 5;
//typedef generic_kernel::kernel<N> kernel_type;
typedef i386_asm_kernel::kernel<N> kernel_type;
typedef symbolic::lint<N,kernel_type> lint;

using namespace std;

int main() {

  int i;
  long long l;
  lint x, y;
  
  cout << " x = ";
  cin >> x;
  cout << " y = ";
  cin >> y;
  cout << " i = ";
  cin >> i;
  cout << " l = ";
  cin >> l;
  cout << " x+y = " << x+y << '\n';
  cout << " x-y = " << x-y << '\n';
  cout << " x*y = " << x*y << '\n';
  cout << " x+i = " << x+i << '\n';
  cout << " x-i = " << x-i << '\n';
  cout << " x*i = " << x*i << '\n';
  cout << " x+l = " << x+l << '\n';
  cout << " x-l = " << x-l << '\n';
  cout << " x*l = " << x*l << '\n';
  cout.flush();

  return 0;

}
