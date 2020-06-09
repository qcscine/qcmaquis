/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2011 by Bela Bauer <bauerb@phys.ethz.ch>
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
 
#ifndef IETL_GMRES_H
#define IETL_GMRES_H
 
#include <vector>
#include <iostream>
 
// Some parts of this code are based on IML++, http://math.nist.gov/iml++/

namespace ietl
{
    namespace detail
    {
        template<class T>
        void GeneratePlaneRotation(T dx, T dy, T & cs, T & sn)
        {
            if (dy == T()) {
                cs = 1; sn = 0;
            } else if (std::abs(dy) > std::abs(dx)) {
                T tmp = dx / dy;
                sn = T(1) / sqrt(T(1)+tmp*tmp);
                cs = tmp*sn;
            } else {
                T tmp = dy / dx;
                cs = T(1) / sqrt(T(1)+tmp*tmp);
                sn = tmp*cs;
            }
        }
        
        template<class T>
        void ApplyPlaneRotation(T & dx, T & dy, T cs, T sn)
        {
            T r0 = cs * dx + sn * dy;
            T r1 = -sn * dx + cs * dy;
            dx = r0;
            dy = r1;
        }
        
        template<class T>
        std::vector<T> Update(boost::numeric::ublas::matrix<T> const & H, std::vector<T> const & S, std::size_t k)
        {
            std::vector<T> y(S.begin(), S.begin()+k);
            for (int i = k-1; i >= 0; --i) {
                y[i] /= H(i,i);
                for (int j = i-1; j >= 0; --j)
                    y[j] -= H(j,i) * y[i];
            }
            return y;
        }
    }
    
    class ietl_gmres
    {
    private:
        std::size_t max_iter;
        bool verbose;
        
    public:
        ietl_gmres(std::size_t max_iter = 100, bool verbose = false)
        : max_iter(max_iter)
        , verbose(verbose) { }
        
        template<class Vector, class Matrix>
        Vector operator()(Matrix const & A,
            Vector const & b,
            Vector const & x0,
            double abs_tol = 1e-6)
        {   
            std::vector<typename Vector::value_type> s(max_iter+1), cs(max_iter+1), sn(max_iter+1);
            std::vector<Vector> v(max_iter+1);
            
            Vector r, w;
            
            mult(A, x0, v[0]);
            r = b - v[0];
            s[0] = two_norm(r);
            
            if (std::abs(s[0]) < abs_tol) {
                if (verbose)
                    std::cout << "Already done with x0." << std::endl;
                return x0;
            }
            
            v[0] = r / s[0];
            
            boost::numeric::ublas::matrix<typename Vector::value_type> H(max_iter+1, max_iter+1);
            std::size_t i = 0;
            
            for ( ; i < max_iter-1; ++i)
            {
                mult(A, v[i], w);
                for (std::size_t k = 0; k <= i; ++k) {
                    H(k,i) = dot(w, v[k]);
                    w -= H(k,i) * v[k];
                }
                
                H(i+1, i) = two_norm(w);
                v[i+1] = w / H(i+1, i);
                
                for (std::size_t k = 0; k < i; ++k)
                    detail::ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
                
                detail::GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
                detail::ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
                detail::ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
                
                if (verbose)
                    std::cout << "GMRES iteration " << i << ", resid = " << std::abs(s[i+1]) << std::endl;
                
                if (std::abs(s[i+1]) < abs_tol)
                    break;
            }
            
            std::vector<typename Vector::value_type> y = detail::Update(H, s, i);
            r = x0;
            for (std::size_t k = 0; k < i; ++k)
                r += y[k] * v[k];
            return r;
        }
    };
}
 
#endif

