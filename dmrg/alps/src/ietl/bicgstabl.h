/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 20XX ??
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


// algorithm according to 
// "Bicgstab(l) For Linear Equations Involving Unsymmetric Matrices With Complex Spectrum (1993)" p. 21 
#ifndef IETL_BICGSTABL_H
#define IETL_BICGSTABL_H
#include <vector>
#include <iostream>
#include <exception>
#include <boost/numeric/ublas/matrix.hpp>
#include <ietl/traits.h>
#include <ietl/complex.h>

namespace ietl {

    namespace ublas = boost::numeric::ublas;
    
    template<class SCALAR>
    class ietl_bicgstabl
    {
    protected:
        unsigned int BICGSTAB_L;
        size_t maxiter;
        bool verbose;
        
    public:
        ietl_bicgstabl(unsigned int b, size_t maxiter = 100, bool verbose = false)
        : BICGSTAB_L(b)
        , maxiter(maxiter)
        , verbose(verbose)
        { }
        
        template <class VECTOR, class MATRIX>
        VECTOR operator() (const MATRIX& A,
                           const VECTOR& b,
                           const VECTOR& x0,
                           typename number_traits<SCALAR>::magnitude_type abs_tol = 1e-6)
        {
            typedef VECTOR vector_t;
            typedef SCALAR scalar_t;
            typedef typename real_type<scalar_t>::type real_t;
            typedef typename number_traits<SCALAR>::magnitude_type magnitude_t;
            typedef std::vector<vector_t> vector_set_t;
            typedef ublas::matrix<scalar_t, ublas::row_major> matrix_t;

            size_t k (0);
            scalar_t rho0(1), rho1, omega(1), alpha(0), beta, gamma;
            vector_t x(x0), rtilde;
            vector_set_t u(BICGSTAB_L+1), r(BICGSTAB_L+1);
            u[0] = vector_t(b.size(),0);

            matrix_t tau(BICGSTAB_L,BICGSTAB_L);
            std::vector<scalar_t> gamma0(BICGSTAB_L), gamma1(BICGSTAB_L), gamma2(BICGSTAB_L), sigma(BICGSTAB_L);

            mult(A, x0, rtilde);
            r[0] = b;
            r[0] -= rtilde;
            rtilde = r[0];

            magnitude_t normr = two_norm(r[0]);

            for(size_t iter = 1; iter <= maxiter && normr > abs_tol; ++iter, k+=BICGSTAB_L)
            {
                rho0 *= -omega;

                for(size_t j = 0; j < BICGSTAB_L; ++j)
                {
                    rho1 = dot(r[j],rtilde);
                    beta = alpha * rho1/rho0;
                    rho0 = rho1;

                    for(size_t i = 0; i < j+1; ++i)
                    {
                        u[i] *= -beta;
                        u[i] += r[i];
                    }

                    mult(A, u[j], u[j+1]);
                    gamma = dot(u[j+1], rtilde);
                    alpha = rho0/gamma;

                    for(size_t i = 0; i < j+1; ++i)
                        r[i] -= alpha * u[i+1];

                    mult(A, r[j], r[j+1]);
                    x += alpha * u[0];
                }

                for(size_t j = 0; j < BICGSTAB_L; ++j)
                {
                    for(size_t i = 0; i < j; ++i)
                    {
                        tau(i,j) = dot(r[j+1],r[i+1]) / sigma[i];
                        r[j+1] -= tau(i,j) * r[i+1];
                    }
                        sigma[j] = dot(r[j+1],r[j+1]);
                        gamma1[j] = dot(r[0],r[j+1]) / sigma[j];
                }

                gamma0[BICGSTAB_L-1] = omega = gamma1[BICGSTAB_L-1];

                for(int j = BICGSTAB_L-2; j >= 0; --j)
                {
                    gamma0[j] = gamma1[j];
                    for(size_t i = j+1; i < BICGSTAB_L; ++i)
                        gamma0[j] -= tau(j,i) * gamma0[i]; //fix?
                }

                for(size_t j = 0; j < BICGSTAB_L-1; ++j)
                {
                    gamma2[j] = gamma0[j+1];
                    for(size_t i = j+1; i < BICGSTAB_L-1; ++i)
                       gamma2[j] += tau(j,i) * gamma0[i+1];
                }

                x += gamma0[0] * r[0];
                r[0] -= gamma1[BICGSTAB_L-1] * r[BICGSTAB_L];
                u[0] -= gamma0[BICGSTAB_L-1] * u[BICGSTAB_L];

                for(size_t j = 0; j < BICGSTAB_L-1; ++j)
                {
                    u[0] -= gamma0[j] * u[j+1];
                    x += gamma2[j] * r[j+1];
                    r[0] -= gamma1[j] * r[j+1];
                }

                normr = two_norm(r[0]);

                if(verbose)
                    std::cout << "BiCGstab("<<BICGSTAB_L<<") Iteration " << iter <<",\t Residual = "<<normr<<"\n";

            }// outer loop

            //if(two_norm(r[0]) > abs_tol)
            //    throw std::runtime_error("BiCGstab failed to converge.");

            return x;
        }
    };

} //end namespace ietl
#endif
