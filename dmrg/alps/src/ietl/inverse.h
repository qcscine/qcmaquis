/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Rene Villiger <rvilliger@smile.ch>,
*                            Prakash Dayal <prakash@comp-phys.org>,
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

/* $Id: inverse.h,v 1.9 2004/02/15 23:30:42 troyer Exp $ */

#ifndef IETL_INVERSE_H
#define IETL_INVERSE_H

#include <ietl/traits.h>
#include <ietl/complex.h>

// REQUIRES class SOLVER with interface
//
//    void operator()(MATRIX matrix, MAGNITUDE_TYPE sigma, VECTOR_TYPE v, VECTOR_TYPE y)
// 
// Solves the equation
//    y = (A - sigma*I)^{-1} * v
//

namespace ietl
{
   template <class MATRIX, class GEN, class SOLVER, class ITER, class VS>
   std::pair<typename vectorspace_traits<VS>::magnitude_type,
             typename vectorspace_traits<VS>::vector_type>
      inverse(const MATRIX& matrix, GEN& gen,const SOLVER& solver, ITER& iter,
                          typename vectorspace_traits<VS>::magnitude_type sigma,
              const VS& vec)
      {
         typedef typename vectorspace_traits<VS>::vector_type vector_type;
         typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
         typedef typename vectorspace_traits<VS>::magnitude_type magnitude_type;
         
         vector_type y   = new_vector(vec);
         vector_type v   = new_vector(vec);
         magnitude_type theta;
         magnitude_type residual;

         // Start with vector y=z, the initial guess
         ietl::generate(y,gen);
         ietl::project(y,vec);
         
         // start iteration loop
         do {
            v=(1./ietl::two_norm(y)) * y;
            try {
                          solver(matrix, sigma, v, y);  // y = (A-\sigma*I)^{-1} v
            }
            catch (...) {
              break; // done with iteration
            }  
                        
            theta = ietl::real(ietl::dot(v,y));
            v = y-theta*v;  // residual = | y - \theta*v |_2
            residual = ietl::two_norm(v);
            ++iter;
            
         // check for convergence
         } while(!iter.finished(residual, theta));

         // accept \lambda = \sigma + 1/\theta  and  x = y/\theta 
         y/=theta;
         theta = sigma + 1./theta;
         
         return std::make_pair(theta, y);
      }
}

#endif
