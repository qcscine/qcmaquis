/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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
 
#ifndef IETL_CG_H
#define IETL_CG_H
 
#include <vector>
#include <iostream>
 
// Note: does not raise any exceptions on failure!
 
namespace ietl
{
    class ietl_cg
    {
    private:
        std::size_t max_iter;
        bool verbose;
        
    public:
        ietl_cg(std::size_t max_iter = 100, bool verbose = false)
        : max_iter(max_iter)
        , verbose(verbose) { }
        
        template<class Vector, class Matrix>
        Vector operator()(Matrix const & A, Vector const & b, Vector const & x0, double abs_tol = 1e-6)
        {   
            std::vector<double> rho(2);
        
            // I assume cheap default-constructibility for now
            Vector x, p, q;
        
            mult(A, x0, x);
            Vector r = b - x;
            
            x = x0;
        
            for (std::size_t k = 1; k < max_iter; ++k)
            {
                rho[0] = rho[1];
                rho[1] = two_norm(r);
                rho[1] *= rho[1];
                if (rho[1] < abs_tol)
                    return x;
                
                if (k == 1)
                    p = r;
                else {
                    double beta = rho[1]/rho[0];
                    p = r + beta*p;
                }
            
                mult(A, p, q);
                double alpha = rho[1]/dot(p, q);
                x += alpha*p;
                r -= alpha*q;
            
                Vector diff = alpha*p;
                double resid = two_norm(diff);
                if (verbose)
                    std::cout << "Iteration " << k << ", resid = " << resid << std::endl;
                if (resid < abs_tol)
                    return x;

            }
            
            return x;
        }
    };
}
 
#endif

