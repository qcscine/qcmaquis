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
 *                            Alberto Baiardi <alberto.baiardi@sns.it>
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
    // +----------------+
    //  Davidson problem
    // +----------------+
    // Class declaration
    template <class MATRIX, class VS>
    class davidson
    {
    public:
        // Interfaces to methods
        typedef typename vectorspace_traits<VS>::vector_type vector_type ;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type ;
        typedef typename std::vector<vector_type> vector_set ;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        davidson(const MATRIX& matrix,
                 const VS& vec);
        davidson(const MATRIX& matrix,
                 const VS& vec,
                 const magnitude_type omega);
        template <class GEN, class PRECOND, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen,
                                                                    PRECOND& mdiag,
                                                                    ITER& iter);
    private:
        // Private attributes
        MATRIX const & matrix_;
        VS vecspace_;
        magnitude_type atol_;
        magnitude_type const omega_ ;
        magnitude_type kappa_ ;
        bool const shift_and_invert_ ;
        // Private methods
        void update_vspace(vector_set& V, vector_set& VA, vector_type& t, std::size_t dim);
    };
    //
    // Constructors
    // ------------
    // NOTE : overloaded to support the shift-and-invert case
    template <class MATRIX, class VS>
    davidson<MATRIX, VS>::davidson(const MATRIX& matrix, const VS& vec) :
            matrix_(matrix),
            vecspace_(vec),
            omega_(0.),
            shift_and_invert_(false),
            kappa_(0.25) {} ;
    template <class MATRIX, class VS>
    davidson<MATRIX, VS>::davidson(const MATRIX& matrix, const VS& vec, const magnitude_type omega ) :
            vecspace_(vec),
            omega_(omega),
            shift_and_invert_(true),
            kappa_(0.25)
    {
        for (int i = 0 ; i < matrix.num_cols ; i++) {
            matrix_(i,i) = omega_ - matrix(i,i);
            for (int j = i+1; j < matrix.num_cols; j++) {
                matrix_(i,j) = -matrix(i,j);
                matrix_(j,i) = matrix_(i,j);
            }
        }
    };
    //
    // Update of the vector space
    // --------------------------
    // TODO better definition of the method depending if shift_and_invert_ is defined or not
    template <class MATRIX, class VS>
    void davidson<MATRIX, VS>::update_vspace(davidson::vector_set &V, davidson::vector_set &VA,
                                             davidson::vector_type &t , std::size_t dim)
    {
        magnitude_type tau = ietl::two_norm(t);
        for (int i = 0; i < dim; i++)
            t -= ietl::dot(V[i], t) * V[i];
        if (ietl::two_norm(t) < this->kappa_ * tau)
            for (int i = 0; i < dim; i++)
                t -= ietl::dot(V[i], t) * V[i];
        V.push_back(t / ietl::two_norm(t));
        VA.resize(V.size());
        ietl::mult(matrix_, V[V.size() - 1], VA[V.size() - 1]);
    }
    //
    // Method used to compute eigenpairs
    // ---------------------------------
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
        std::size_t iter_dim = 0 ;
        unsigned int i,j;
        magnitude_type theta, tau;
        magnitude_type rel_tol;
        atol_ = iter.absolute_tolerance();
        // Start with t=v_o, starting guess
        ietl::generate(t,gen);
        ietl::project(t,vecspace_);
        // Start iteration
        do
        {
            // Modified Gram-Schmidt Orthogonalization with Refinement
            update_vspace(V, VA, t, iter_dim);
            iter_dim = V.size() ;
            matrix_t M(iter_dim, iter_dim), Mevecs(iter_dim, iter_dim);
            std::vector<magnitude_type> Mevals(iter_dim);
            for (i = 0; i < iter_dim; ++i)
                for (j = i; j < iter_dim; ++j)  {
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
            if (iter.finished(ietl::two_norm(r), Mevals[0]))
                break;
            mdiag.precondition(r, V[0], Mevals[0]);
            std::swap(t,r);
            t /= ietl::two_norm(t);
            //if (V.size() >= 20)
            //{
            //    V.resize(2);
            //    VA.resize(2);
            //}
        } while (true);
        // accept lambda=theta and x=u
        return std::make_pair(theta, u);
    }
    //
    // +-------------------------+
    //  Modified Davidson problem
    // +-------------------------+
    // Class declaration
    template <class MATRIX, class VS, class SymmGroup, class OtherMatrix>
    class davidson_modified
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        typedef typename std::vector<vector_type> vector_set ;
        typedef typename partial_overlap<OtherMatrix, SymmGroup>::partial_overlap partial_overlap ;
        davidson_modified(const MATRIX& matrix,
                          const VS& vec,
                          const double& omega,
                          const partial_overlap& poverlap,
                          const int& site);
        template <class GEN, class PRECOND, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen,
                                                                    PRECOND& mdiag,
                                                                    ITER& iter);
    private:
        MATRIX const & matrix_;
        VS vecspace_;
        magnitude_type atol_;
        double omega_ ;
        void update_vspace(vector_set& V,  vector_set& VA, vector_set& V2, vector_type& t , std::size_t dim);
        vector_type apply_operator(const vector_type& x);
        partial_overlap poverlap_ ;
        int site_ ;
    };
    // Constructor
    template <class MATRIX, class VS, class SymmGroup, class OtherMatrix>
    davidson_modified<MATRIX, VS, SymmGroup, OtherMatrix>::davidson_modified(const MATRIX& matrix,
                                                                             const VS& vec,
                                                                             const double& omega,
                                                                             const partial_overlap& poverlap,
                                                                             const int& site) :
            vecspace_(vec),
            matrix_(matrix),
            omega_(omega),
            poverlap_(poverlap),
            site_(site)
    { }
    template <class MATRIX, class VS, class SymmGroup, class OtherMatrix>
    typename davidson_modified<MATRIX, VS, SymmGroup, OtherMatrix>::vector_type davidson_modified<MATRIX, VS, SymmGroup, OtherMatrix>::apply_operator(const davidson_modified::vector_type& x)
    {
        vector_type tmp ;
        ietl::mult(matrix_, x, tmp);
        tmp *= -1. ;
        tmp += omega_*x ;
        return tmp ;
    };
    template <class MATRIX, class VS, class SymmGroup, class OtherMatrix>
    void davidson_modified<MATRIX, VS, SymmGroup, OtherMatrix>::update_vspace(davidson_modified::vector_set& V,  davidson_modified::vector_set& VA,
                                                                              davidson_modified::vector_set& V2, davidson_modified::vector_type& t ,
                                                                              std::size_t dim)
    {
        magnitude_type tau = ietl::two_norm(t);
        vector_type tA ;
        for (int i = 0; i < dim; i++)
            t -= ietl::dot(V[i], t) * V[i];
        t /= ietl::two_norm(t);
        V.push_back(t);
        tA = apply_operator(t) ;
        for (int i = 0; i < dim; i++) {
            t -= ietl::dot(VA[i], tA) * V2[i];
            tA -= ietl::dot(VA[i], tA) * VA[i];
        }
        V2.push_back(t/ietl::two_norm(tA));
        VA.push_back(tA/ietl::two_norm(tA));
    }
    // Method to actually diagonalize the matrix
    template <class MATRIX, class VS, class SymmGroup, class OtherMatrix>
    template <class GEN, class PRECOND, class ITER>
    std::pair<typename davidson_modified<MATRIX, VS, SymmGroup, OtherMatrix>::magnitude_type,
              typename davidson_modified<MATRIX, VS, SymmGroup, OtherMatrix>::vector_type>
    davidson_modified<MATRIX, VS, SymmGroup, OtherMatrix>::calculate_eigenvalue(const GEN& gen,
                                                                   PRECOND& mdiag,
                                                                   ITER& iter)
    {
        // Type definition
        typedef alps::numeric::matrix<scalar_type> matrix_t;
        vector_type t  = new_vector(vecspace_);
        vector_type tA = new_vector(vecspace_);
        vector_type tB = new_vector(vecspace_);
        vector_type u  = new_vector(vecspace_);
        vector_type uA = new_vector(vecspace_);
        vector_type uB = new_vector(vecspace_);
        vector_type r, r2  = new_vector(vecspace_);
        std::vector<scalar_type> s(iter.max_iterations());
        std::vector<vector_type> V;
        std::vector<vector_type> V2;
        std::vector<vector_type> VA;
        unsigned int i,j;
        magnitude_type theta, tau;
        magnitude_type kappa = 0.25;
        magnitude_type rel_tol;
        atol_ = iter.absolute_tolerance();
        std::size_t iter_dim = 0 ;
        // Start with t=v_o, starting guess
        ietl::generate(t,gen);
        ietl::project(t,vecspace_);
        // Magnitude type
        magnitude_type shift = omega_ ;
        // Start iteration
        do {
            update_vspace(V, VA, V2, t, iter_dim);
            iter_dim = V2.size();
            matrix_t M(iter_dim, iter_dim), Mevecs(iter_dim, iter_dim);
            std::vector<magnitude_type> Mevals(iter_dim);
            for (i = 0; i < iter_dim; ++i)
                for (j = i; j < iter_dim; ++j)
                {
                    M(i,j) = ietl::dot(V2[i], VA[j]);
                    M(j,i) = M(i,j);
                }
            boost::numeric::bindings::lapack::heevd('V', M, Mevals);
            Mevecs = M ;
            u = V2[0]*Mevecs(0,0) ;
            for (i = 1; i < iter_dim; ++i)
                u += V2[i] * Mevecs(i, 0);
            uA = VA[0]*Mevecs(0,0);
            for (i = 1; i < iter_dim; ++i)
                uA += VA[i] * Mevecs(i,0);
            uA /= ietl::two_norm(u) ;
            u  /= ietl::two_norm(u) ;
            magnitude_type energy = ietl::dot(u,uA)  ;
            r2 = uA - energy*u ;
            theta = energy ;
            // if (|r|_2 < \epsilon) stop
            ++iter;
            if (iter.finished(ietl::two_norm(r2), theta))
                break;
            mdiag.precondition(r2, u, theta);
            std::swap(t,r2);
            //if (V.size() >= 20)
            //{
            //    V.resize(2);
            //    V2.resize(2);
            //    VA.resize(2);
            //}
        } while (true);
    // accept lambda=theta and x=u
    return std::make_pair(shift-theta, u);
    }
}
#endif
