/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id: blitz.h,v 1.10 2003/09/05 08:12:38 troyer Exp $ */

#ifndef IETL_INTERFACE_BLITZ_H
#define IETL_INTERFACE_BLITZ_H

#include <ietl/complex.h>
#include <blitz/array.h>
#include <blitz/vecglobs.h>
#include <ietl/traits.h>


namespace ietl {
  
  template < class Cont>
    void clear(Cont& c) {
    std::fill(c.begin(),c.end(),0.);
  }

  template < class Cont, class Gen> 
    void generate(Cont& c, const Gen& gen) {
    std::generate(c.begin(),c.end(),gen);
  }

  template <class T, int D>
  typename number_traits<T>::magnitude_type two_norm(const blitz::Array<T,D>& v) {
    return std::sqrt(ietl::real(dot(v,v)));
  }

  template <class T, int D>
  T dot(const blitz::Array<T,D>& x, const blitz::Array<T,D>& y) {
    return blitz::sum(x*y);
  }

  template <class T, int D>
  T dot(const blitz::Array<std::complex<T>,D>& x, const blitz::Array<std::complex<T>,D>& y) {
    
    return blitz::sum(blitz::conj(x)*y);
  }

  template <class T, int D>
  void copy(const blitz::Array<T,D>& x, blitz::Array<T,D>& y) {
    y=x;
    y.makeUnique();
  }
}

namespace std {  
  template <class T, int D>
  void swap(blitz::Array<T,D>& x, blitz::Array<T,D>& y) {
    blitz::cycleArrays(x,y);
  }

}

#endif // IETL_INTERFACE_BLITZ_H

