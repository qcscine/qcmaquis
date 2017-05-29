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
    //
    // CLASS DAVIDSON
    // --------------
    //
    // Has private attributes:
    // - MATRIX : the matrix to be diagonalized
    // - VS : a vector space
    // - atol: tolerance for the Davidson algorithm
    // - DesiredEigenvalue : which eigenvalue to compute
    //
    template <class MATRIX, class VS>
    class davidson
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        
        
        davidson(const MATRIX& matrix, 
                 const VS& vec,
                 DesiredEigenvalue desired = Largest);
        
        template <class GEN, class SOLVER, class PRECOND, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen, 
                                                                    SOLVER& solver,
                                                                    PRECOND& mdiag,
                                                                    ITER& iter);
    private:
        MATRIX const & matrix_;
        VS vecspace_;
        magnitude_type atol_;
        DesiredEigenvalue desired_;
    };
    
    template <class MATRIX, class VS>
    davidson<MATRIX, VS>::davidson(const MATRIX& matrix, const VS& vec, DesiredEigenvalue desired) : 
    matrix_(matrix),
    vecspace_(vec),
    desired_(desired)
    {}
    
    template <class MATRIX, class VS> 
    template <class GEN, class SOLVER, class PRECOND, class ITER>
    //
    // BRIEF DESCRIPTION OF THE DAVIDSON ALGORITHM
    // The Davidson algorithm is used to compute eigenvectors or large matrices by expanding the
    // eigenvector in a (relatively) small vector space. The idea:
    // 1 - start from a guess psi_0
    // 2 - compute the error (H-E)psi_0, where E is the expectation value of the Hamiltonian over
    //     psi_0
    // 3 - multiply the error by (diag(H)-E)^{-1} (just like a preconditioner)
    // 4 - take psi_1 as the resulting vector, orthogonalize (here a refined Gram-Schmidt algorithm
    //     is employed) wrt to psi_0
    // 5 - restart from 2 and increase the subspace until convergence
    //
    std::pair<typename davidson<MATRIX,VS>::magnitude_type, typename davidson<MATRIX, VS>::vector_type> 
    davidson<MATRIX, VS>::calculate_eigenvalue(const GEN& gen,
                                               SOLVER& solver,
                                               PRECOND& mdiag,
                                               ITER& iter)
    {
        typedef alps::numeric::matrix<scalar_type> matrix_t;
        // t is the trial vector, updated at each cycle
        vector_type t  = new_vector(vecspace_);
        //ALB
        vector_type tA = new_vector(vecspace_);
        vector_type u  = new_vector(vecspace_);
        vector_type uA = new_vector(vecspace_);
        vector_type r  = new_vector(vecspace_);
        // V collects the vectors at the previous cycles, VA collects the
        // results of the application of the Hamiltonian matrix to V
        // TODO check if s is actually used
        std::vector<scalar_type> s(iter.max_iterations());
        std::vector<vector_type> V;
        std::vector<vector_type> VA;

        unsigned int i,j;
        magnitude_type theta, tau;
        magnitude_type kappa = 0.25;
        atol_ = iter.absolute_tolerance();
        // -- Initialization --
        ietl::generate(t,gen);
        ietl::project(t,vecspace_);
        // -- Main iteration --
        do
        {
            // Modified Gram-Schmidt Orthogonalization with Refinement
            // NOTE: V contains the vectors determined in the previous cycles
            magnitude_type scal ;
            tau = ietl::two_norm(t);
            ietl::mult( matrix_ , t , tA ) ;
            for (i = 0; i < VA.size(); i++) {
                scal = ietl::dot(30.*t - tA, 30.*V[i] - VA[i]);
                tA -= scal * (30.*V[i] - VA[i]) ;
                t  -= scal * V[i];
            }
            // TODO ALB check later how to generalize this part
            //if (ietl::two_norm(t) < kappa * tau)
            //    for (i = 0; i < V.size(); i++)
            //        t -= ietl::dot(V[i], t) * V[i];
            // Project out orthogonal subspace
            //ietl::project(t, vecspace_);
            // Update of V and VA with the new vector
            magnitude_type normalize ;
            normalize = ietl::two_norm(30.*t - tA) ;
            V.push_back(t/normalize);
            VA.resize(V.size());
            // Put in VA the result of the multiplication of V times matrix_
            // (only the last vector is updated)
            ietl::mult(matrix_, V[V.size() - 1], VA[V.size() - 1]);
            // Variable declaration
            std::size_t iter_dim = V.size();
            matrix_t M(iter_dim, iter_dim), Mevecs(iter_dim, iter_dim);
            std::vector<magnitude_type>     Mevals(iter_dim);
            // Main cycle: compute the representation of M in the V space
            for (i = 0; i < iter_dim; ++i)
                for (j = i; j < iter_dim; ++j) {
                    M(i,j) = 30.*ietl::dot(V[i],V[j]) - ietl::dot(V[i], VA[j]) ;
                    M(j,i) = M(i,j);
                }
            // Diagonalize (NOTE: the matrix is destroyed in output)
            // TODO check if stuff is sorted in output
            boost::numeric::bindings::lapack::heevd('V', M, Mevals);
            Mevecs = M;
            // Bring back to the original space
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
            // Compute the error for the lowest-energy solution
            // VA = H*psi and V*MEVals = E*psi
            r = (30.*V[0] - VA[0]) - V[0] * Mevals[0] ;
            theta = Mevals[0] ;
            u = V[0];
            // if (|r|_2 < \epsilon) stop
            ++iter;
            if (iter.finished(ietl::two_norm(r), Mevals[0]))
                break;
            // Take the remainder and multiply it by (diag(H)-E) and uses it as
            // a new vector. Iterate until convergence
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
        return std::make_pair(30.-1.0/theta, u);
    }
}
#endif
