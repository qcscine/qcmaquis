/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef IETL_BCONJGRAD_STAB_ALB_H
#define IETL_BCONJGRAD_STAB_ALB_H

#include <vector>
#include <iostream>

//
// This code is based on IML++, http://math.nist.gov/iml++/
// Modified by Alberto Baiardi to fix an inconsistency in the
// definition of the indexes in the original alps version
//

namespace ietl
{
    //
    // GMRES_GENERAL CLASS
    // -------------------
    // General class with the GMRES solver. The definition of the matrix of the linear system
    // is left free, since it differs for standard and shift-and-invert problem
    template<class Matrix, class Vector, class VectorSpace>
    class biconjgrad_stab_general {
        //
    public:
        // Types definition
        typedef typename Vector::value_type                            scalar_type ;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        typedef typename std::size_t                                   size_t ;
        typedef typename boost::numeric::ublas::matrix< scalar_type >  matrix_scalar ;
        typedef typename std::vector<scalar_type>                      vector_scalar ;
        typedef typename std::vector<Vector>                           vector_space ;
        typedef typename std::pair<Vector , Vector >                   pair_vectors ;
        typedef typename std::vector< std::vector<Vector> >            vector_ortho_vec ;
        //
        // Public methods
        // --------------
        // Constructor
        biconjgrad_stab_general(Matrix const & A,
                                Vector const & u,
                                VectorSpace const & vs,
                                magnitude_type const & theta,
                                vector_ortho_vec const & ortho_vec,
                                size_t const& n_root,
                                std::size_t max_iter = 100,
                                bool verbose = false)
                : max_iter(max_iter), verbose(verbose), A_(A), u_(u), theta_(theta), vs_(vs), n_root_(n_root), ortho_vec_(ortho_vec)
        { }
        // () operator. Returns the solution of the GMRES problem
        Vector operator()(Vector const & b,
                          Vector const & x0,
                          double abs_tol = 1e-6,
                          double rel_tol = 1e-3)
        {
            // Types definiton
            Vector x = x0 ;     //  Approximation of the solution at the current iteration
            Vector p , Ap ;     //  Conjugate vector at the current iteration, and the results of the application of A to p
            Vector s , As ;
            Vector h ;
            Vector r, r0 ;      //  Error vectors
            Vector jnk ;        //  Scratch vector
            scalar_type alpha ; //  Coefficient to use to update x
            scalar_type beta  ;
            scalar_type omega ;
            scalar_type normb = ietl::two_norm(b) , normr ;
            scalar_type rho_pr, rho_act ;
            // Initialization (0-th iteration)
            Vector init = x0 ;
            jnk = apply(init) ;
            r0 = b - jnk ;
            r  = r0 ;
            p  = 0.*r ;
            Ap = 0.*r ;
            rho_pr = 1. ;
            alpha  = 1. ;
            omega  = 1. ;
            for (std::size_t i = 0 ; i < max_iter ; ++i)
            {
                rho_act = ietl::dot(r0, r) ;
                beta    = (rho_act*alpha)/(rho_pr*omega) ;
                jnk     = beta*p ;
                p       = r + jnk - beta*omega*Ap ;
                Ap      = apply(p) ;
                alpha   = rho_act / ietl::dot(r0, Ap) ;
                h       = x + alpha*p ;
                jnk     = b - apply(h) ;
                if (ietl::two_norm(jnk) < abs_tol)
                    return h ;
                s  = r - alpha*Ap ;
                As = apply(s) ;
                omega = ietl::dot(As,s) / ietl::dot(As,As) ;
                x     = h + omega*s ;
                jnk   = b - apply(x) ;
                if (verbose)
                    std::cout << "Conj. Grad. iteration " << i << ", Abs. Err. = " << ietl::two_norm(jnk)
                                                               << ", Rel. Err. = " << ietl::two_norm(jnk)/normb
                                                               << std::endl ;
                if ( maquis::real(ietl::two_norm(jnk)/normb) < rel_tol || maquis::real(ietl::two_norm(jnk)) < abs_tol ) {
                    return x;
                } else {
                    rho_pr = rho_act;
                    r = s - omega * As;
                }
            }
            return x ;
        }
    protected:
        // Compute the result of the application of the operator
        // This function is left virtual since its definition is different for standard
        // and shift-and-invert problems
        virtual Vector apply(Vector& input) {} ;
        // Orthogonalization routine
        void orthogonalize_simple(Vector& input)
        {
            for (typename vector_ortho_vec::iterator it = ortho_vec_.begin(); it != ortho_vec_.end(); it++)
                if (maquis::real(ietl::dot((*it)[0], (*it)[0])) > 1.0E-15)
                    input -= ietl::dot((*it)[0], input) * (*it)[0] / ietl::dot((*it)[0], (*it)[0]) ;
            return ;
        }
        //
        //
        // Private attributes
        // ------------------
        // A_       --> Matrix to be diagonalized
        // u_       --> Vector used to determine the orthogonal subspace
        // max_iter --> Maximum number of iterations for the GMRES algorithm
        // vebose   --> Controls the output printing
        bool verbose ;
        double theta_ ;
        Matrix A_ ;
        std::size_t max_iter , n_root_ ;
        Vector u_ ;
        VectorSpace vs_ ;
        vector_ortho_vec ortho_vec_ ;
    };
    //
    // CONJGRAD_STANDARD OBJECT
    // ---------------------
    // Conjugate gradient solver object for the standard problem
    template<class Matrix, class Vector, class VectorSpace>
    class biconjgrad_stab_standard : private biconjgrad_stab_general<Matrix, Vector, VectorSpace>
    {
    public:
        typedef typename Vector::value_type                scalar_type ;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        typedef biconjgrad_stab_general<Matrix, Vector, VectorSpace> base ;
        typedef typename base::size_t                      size_t ;
        typedef typename base::vector_ortho_vec            vector_ortho_vec ;
        //
        biconjgrad_stab_standard(Matrix const & A,
                                 Vector const & u,
                                 VectorSpace const & vs,
                                 magnitude_type const & theta,
                                 const vector_ortho_vec& ortho_vec,
                                 size_t const& n_root,
                                 size_t max_iter,
                                 bool verbose)
                : base::biconjgrad_stab_general(A, u, vs, theta, ortho_vec, n_root, max_iter, verbose) { }
        // Private attributes
        using base::A_ ;
        using base::n_root_ ;
        using base::ortho_vec_ ;
        using base::theta_ ;
        using base::u_ ;
        using base::vs_ ;
        using base::operator() ;
    private:
        Vector apply(Vector& input){
            // Initialization
            Vector t, t2, t3, y;
            // t2 = (1-uu*) x
            double ust = maquis::real(dot(u_, input) / ietl::dot(u_, u_)) ;
            t2 = input - ust * u_ ;
            this->orthogonalize_simple(t2) ;
            // y = (A-theta*1) t2
            ietl::mult(A_, t2, t3, n_root_);
            y = t3 - theta_ * t2;
            this->orthogonalize_simple(y) ;
            // t = (1-uu*) y
            ust = maquis::real(dot(u_, y) / ietl::dot(u_, u_)) ;
            t = y - ust * u_ ;
            y = t ;
            // Finalization
            return y ;
        }
    };
}

#endif
