/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2015 by Rene Villiger <rvilliger@smile.ch>,
 *                            Prakash Dayal <prakash@comp-phys.org>,
 *                            Matthias Troyer <troyer@comp-phys.org>
 *                            Bela Bauer <bauerb@phys.ethz.ch>
 *                            Sebastian Keller <sebkelle@phys.ethz.ch>
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

// Davidson algorithm adapted from the IETL Jacobi-Davidson implementation in ALPS

#ifndef IETL_DAVIDSON_H
#define IETL_DAVIDSON_H

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h> 
 
#include <ietl/cg.h>
#include <ietl/gmres.h>

#include <complex>
#include <vector>

#include <boost/function.hpp>

namespace ietl
{
    // +-------------------------+
    //  Standard Davidson problem
    // +-------------------------+
    // Class declaration
    template <class MATRIX, class VS>
    class davidson
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        davidson(const MATRIX& matrix,
                 const VS& vec);
        template <class GEN, class PRECOND, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen,
                                                                    PRECOND& mdiag,
                                                                    ITER& iter);
    private:
        MATRIX const & matrix_;
        VS vecspace_;
        magnitude_type atol_;
    };
    // Class
    template <class MATRIX, class VS>
    davidson<MATRIX, VS>::davidson(const MATRIX& matrix, const VS& vec) :
            matrix_(matrix),
            vecspace_(vec)
    {}
    // Diagonalization routine
    template <class MATRIX, class VS>
    template <class GEN, class PRECOND, class ITER>
    std::pair<typename davidson<MATRIX,VS>::magnitude_type, typename davidson<MATRIX, VS>::vector_type>
    davidson<MATRIX, VS>::calculate_eigenvalue(const GEN& gen,
                                               PRECOND& mdiag,
                                               ITER& iter)
    {
        // Initialization
        typedef alps::numeric::matrix<scalar_type> matrix_t;
        vector_type t  = new_vector(vecspace_);
        vector_type u  = new_vector(vecspace_);
        vector_type r  = new_vector(vecspace_);
        std::vector<scalar_type> s(iter.max_iterations());
        std::vector<vector_type> V;
        std::vector<vector_type> VA;
        unsigned int i,j;
        magnitude_type theta, tau;
        magnitude_type kappa = 0.25;
        magnitude_type rel_tol;
        atol_ = iter.absolute_tolerance();
        // Start with t=v_o, starting guess
        ietl::generate(t,gen);
        ietl::project(t,vecspace_);
        // Start iteration
        do
        {
            // Modified Gram-Schmidt Orthogonalization with Refinement
            tau = ietl::two_norm(t);
            for (i = 0; i < V.size(); i++)
                t -= ietl::dot(V[i], t) * V[i];
            if (ietl::two_norm(t) < kappa * tau)
                for (i = 0; i < V.size(); i++)
                    t -= ietl::dot(V[i], t) * V[i];
            // Project out orthogonal subspace
            ietl::project(t, vecspace_);
            // v_m = t / |t|_2,  v_m^A = A v_m
            V.push_back(t/ietl::two_norm(t));
            VA.resize(V.size());
            ietl::mult(matrix_, V[V.size() - 1], VA[V.size() - 1]);
            std::size_t iter_dim = V.size();
            matrix_t M(iter_dim, iter_dim), Mevecs(iter_dim, iter_dim);
            std::vector<magnitude_type> Mevals(iter_dim);
            for (i = 0; i < iter_dim; ++i)
                for (j = i; j < iter_dim; ++j)
                {
                    M(i,j) = ietl::dot(V[i], VA[j]);
                    M(j,i) = M(i,j);
                }
            boost::numeric::bindings::lapack::heevd('V', M, Mevals);
            Mevecs = M;
            std::vector<vector_type> Vp = V, VAp = VA;
            for (i = 0; i < iter_dim; ++i)
            {
                V[i] *= Mevecs(i,i);
                VA[i] *= Mevecs(i,i);
            }
            for (i = 0; i < iter_dim; ++i)
                for (j = 0; j < iter_dim; ++j)
                    if(i != j)
                    {
                        V[j] += Vp[i] * Mevecs(i,j);
                        VA[j] += VAp[i] * Mevecs(i,j);
                    }
            r = VA[0] - V[0] * Mevals[0];
            theta = Mevals[0];
            u = V[0];
            // if (|r|_2 < \epsilon) stop
            ++iter;
            std::cout << two_norm(r)/Mevals[0] << std::endl ;
            if (iter.finished(ietl::two_norm(r), Mevals[0]))
                break;
            mdiag.precondition(r, V[0], Mevals[0]);
            std::swap(t,r);
            t /= ietl::two_norm(t);
            if (V.size() >= 20)
            {
                V.resize(2);
                VA.resize(2);
            }
        } while (true);
        // accept lambda=theta and x=u
        return std::make_pair(theta, u);
    }
    //
    // +-------------------------+
    //  Modified Davidson problem
    // +-------------------------+
    // Class declaration
    template <class MATRIX, class VS>
    class davidson_modified
    {
    public:
        // TODO generalize the setting of the type of omega
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        davidson_modified(const MATRIX& matrix,
                          const VS& vec,
                          double omega);
        template <class GEN, class PRECOND, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen,
                                                                    PRECOND& mdiag,
                                                                    ITER& iter);
    private:
        MATRIX const & matrix_;
        VS vecspace_;
        magnitude_type atol_;
        double omega_ ;
    };
    // Constructor
    template <class MATRIX, class VS>
    davidson_modified<MATRIX, VS>::davidson_modified(const MATRIX& matrix, const VS& vec, const double omega) :
            matrix_(matrix),
            vecspace_(vec),
            omega_(omega)
    {}
    // Method to actually diagonalize the matrix
    template <class MATRIX, class VS>
    template <class GEN, class PRECOND, class ITER>
    std::pair<typename davidson_modified<MATRIX,VS>::magnitude_type, typename davidson_modified<MATRIX, VS>::vector_type>
    davidson_modified<MATRIX, VS>::calculate_eigenvalue(const GEN& gen,
                                                        PRECOND& mdiag,
                                                        ITER& iter)
    {
        // Type definition
        typedef alps::numeric::matrix<scalar_type> matrix_t;
        vector_type t  = new_vector(vecspace_);
        vector_type tA = new_vector(vecspace_);
        vector_type u  = new_vector(vecspace_);
        vector_type uA = new_vector(vecspace_);
        vector_type uB = new_vector(vecspace_);
        vector_type r  = new_vector(vecspace_);
        std::vector<scalar_type> s(iter.max_iterations());
        std::vector<vector_type> V;
        std::vector<vector_type> VA;
        unsigned int i,j;
        magnitude_type theta, tau;
        magnitude_type kappa = 0.25;
        magnitude_type rel_tol;
        atol_ = iter.absolute_tolerance();
        // Start with t=v_o, starting guess
        ietl::generate(t,gen);
        ietl::project(t,vecspace_);
        // Magnitude type
        magnitude_type shift = omega_ ;
        // Start iteration
        do {
            // Modified Gram-Schmidt Orthogonalization with Refinement
            tau = ietl::two_norm(t);
            ietl::mult(matrix_ , t , tA);
            tA *= -1 ;
            tA += shift*t ;
            for (i = 0; i < VA.size(); i++) {
                t -= ietl::dot(VA[i], tA) * V[i];
                tA -= ietl::dot(VA[i], tA) * VA[i];
            }
            //if (ietl::two_norm(t) < kappa * tau)
            //    for (i = 0; i < V.size(); i++)
            //        t -= ietl::dot(V[i], t) * V[i];
            // Project out orthogonal subspace
            //ietl::project(t, vecspace_);
            // v_m = t / |t|_2,  v_m^A = A v_m
            V.push_back(t/ietl::two_norm(tA));
            VA.push_back(tA/ietl::two_norm(tA));
            std::size_t iter_dim = V.size();
            matrix_t M(iter_dim, iter_dim), Mevecs(iter_dim, iter_dim);
            std::vector<magnitude_type> Mevals(iter_dim);
            for (i = 0; i < iter_dim; ++i)
                for (j = i; j < iter_dim; ++j)
                {
                    M(i,j) = ietl::dot(V[i], VA[j]);
                    M(j,i) = M(i,j);
                }
            boost::numeric::bindings::lapack::heevd('V', M, Mevals);
            Mevecs = M;
            std::vector<vector_type> Vp = V, VAp = VA;
            for (i = 0; i < iter_dim; ++i)
            {
                V[i] *= Mevecs(i,i);
                VA[i] *= Mevecs(i,i);
            }
            for (i = 0; i < iter_dim; ++i)
                for (j = 0; j < iter_dim; ++j)
                    if(i != j)
                    {
                        V[j] += Vp[i] * Mevecs(i,j);
                        VA[j] += VAp[i] * Mevecs(i,j);
                    }
            u = V[0]/ietl::two_norm(V[0]);
            ietl::mult(matrix_, u, uA);
            uB = shift*u - uA;
            r = ( uB - u/Mevals[0] );
            theta = 1./Mevals[0];
            // if (|r|_2 < \epsilon) stop
            ++iter;
            std::cout << ietl::two_norm(r) << std::endl ;
            if (iter.finished(ietl::two_norm(r), 1./Mevals[0]))
                break;
            mdiag.precondition(r, u, 1./Mevals[0]);
            std::swap(t,r);
            t /= ietl::two_norm(t);
            if (V.size() >= 20)
            {
                V.resize(2);
                VA.resize(2);
            }
        } while (true);
    // accept lambda=theta and x=u
    return std::make_pair(shift-theta, u);
    }
}
#endif
