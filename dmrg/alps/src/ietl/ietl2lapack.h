/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2010 by Prakash Dayal <prakash@comp-phys.org>,
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

/* $Id: ietl2lapack.h,v 1.9 2003/12/05 09:24:01 tprosser Exp $ */

#ifndef IETL_LAPACK_H
#define IETL_LAPACK_H

#include <boost/numeric/bindings/lapack/driver/stev.hpp>
#include <boost/numeric/bindings/lapack/computational/steqr.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/heev.hpp>
#include <complex>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <cassert>

#undef minor

namespace ietl {
  template <class T>
    T* get_data(const std::vector<T>& v) {
    if (v.empty())
      return 0;
    else
      return const_cast<T*>(&v[0]);
  }
}
 
namespace ietl_lapack_dispatch {  

  inline void stev(const char& jobz, fortran_int_t n, double dd[],double de[],
                   double dz[],fortran_int_t ldz,fortran_int_t info) {
    double* dwork = new double[2*n -2];
    LAPACK_DSTEV(&jobz,&n,dd,de,dz,&ldz,dwork,&info);
    delete[] dwork;
    if (info)
      throw std::runtime_error("Error return from dstev");
  }


  inline void stev(const char& jobz, fortran_int_t n, double dd[],double de[],
                   std::complex<double> dz[],fortran_int_t ldz,fortran_int_t info) {
    double* dwork = new double[2*n -2];
    LAPACK_ZSTEQR(&jobz,&n,dd,de,dz,&ldz,dwork,&info);
    delete[] dwork;
    if (info)
      throw std::runtime_error("Error return from zsteqr");
  }
  
} // ietl2lapack_dispatch ends here.

namespace ietl2lapack {
      
  template<class Vector>
    fortran_int_t stev(const Vector& alpha, const Vector& beta, Vector& eval, fortran_int_t n) {  
    if (n==0) n = alpha.size();
    std::copy(alpha.begin(),alpha.begin() + n, eval.begin()); 
    assert(eval.size() >= n);
    assert(alpha.size() >= n);
    assert(beta.size() >= n);
    Vector beta_tmp(n);
    std::copy(beta.begin(),beta.begin() + n, beta_tmp.begin()); 
    Vector z; // not referenced
    char _jobz = 'N';
    fortran_int_t _info=0;    
    fortran_int_t _ldz = 1; 
    ietl_lapack_dispatch::stev(_jobz, n, ietl::get_data(eval), ietl::get_data(beta_tmp), ietl::get_data(z),_ldz, _info);
    return _info;
  }
  
  template<class Vector, class FortranMatrix>
    int stev(const Vector& alpha, const Vector& beta, Vector& eval, FortranMatrix& z, fortran_int_t n) {  
    if (n==0) n = alpha.size();
    std::copy(alpha.begin(),alpha.begin() + n, eval.begin()); 
        assert(eval.size()>=static_cast<std::size_t>(n));
    assert(alpha.size()>=static_cast<std::size_t>(n));
    assert(beta.size()>=static_cast<std::size_t>(n)); 
    Vector beta_tmp(n);
    std::copy(beta.begin(),beta.begin() + n, beta_tmp.begin()); 
    char _jobz;
    fortran_int_t _info=0;
    _jobz = 'V';
    fortran_int_t _ldz = z.minor(); 
    ietl_lapack_dispatch::stev(_jobz, n, ietl::get_data(eval), ietl::get_data(beta_tmp),z.data(),_ldz, _info);
    return _info;
  }
  
  
} 
#endif
