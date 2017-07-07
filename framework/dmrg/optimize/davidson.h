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

// +--------------+
//  DAVIDSON CLASS
// +--------------+
// This is a general class for Davidson-type eigensolver.
// The templates arguments are MATRIX and VS, that are usually:
// MATRIX    : a SiteProblem object (see optimize.h for additional details)
// VS        : a VectorSpace object, including a MPSTensor and several other vectors
//             (for excited states orthogonalization)
//
// Includes the following attributes, that are common to all the Davidson eigensolvers;
// matrix_   : the matrix representation of the operator for the site where the optimization is
//             carried out
// vecspace_ : the vector space where the optimization is carried out
// atol_     : tolerance for assessing convergence of the optimization
// Other attributes, that are specific of other type of eigensolvers (such as a shift omega or a
// state to target) are defined in the inherited classes.
//
// The methods are the following:
// contructor           : standard constructor
// calculate_eigenvalue : method to compute eigenvalue, uses virtual functions
// update_vspace        : virtual protected function, defined when virtual class is inherited
// apply_operator       : virtual protected function, defined when virtual class is inherited
// precondition         : virtual protected function, for the guess iteration of each Davidson step

namespace ietl
{
    template <class MATRIX, class VS>
    class davidson
    {
    public:
        // Types definition
        typedef typename vectorspace_traits<VS>::vector_type 		      vector_type ;
        typedef typename vectorspace_traits<VS>::scalar_type 		      scalar_type ;
        typedef typename vector_type::bm_type 		      		          bm_type ;
        typedef typename std::vector<vector_type> 	     		          vector_set ;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        typedef typename alps::numeric::matrix<scalar_type> 		      matrix_numeric ;
        typedef typename std::size_t 					                  size_t ;
        // Constructor and destructor
        davidson(const MATRIX& matrix, const VS& vec, const int& site, const int& nmin, const int& nmax);
        virtual ~davidson() {};
        // Public method to compute eigenvalue
        template <class GEN, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen,
                                                                    ITER& iter);
    protected:
        // Virtual Methods, defined in the derived classes
        virtual magnitude_type return_final(const magnitude_type& eigval) {} ;
        virtual vector_type apply_operator(const vector_type& x) {} ;
        virtual vector_type finalize_iteration(const vector_type& u, const vector_type& r, const size_t& n_restart,
                                               size_t& iter_dim, vector_set& v2, vector_set& VA) {} ;
        virtual void precondition(vector_type& r, const vector_type& V, const magnitude_type theta ) {} ;
	    virtual void select_eigenpair(const vector_set& V, const vector_set& VA, const matrix_numeric& eigvecs,
		                 		      const size_t& i, vector_type& u, vector_type& uA) {} ;
        virtual void update_vspace(vector_set& V, vector_set& VA, vector_type& t, std::size_t dim) {} ;
        // Attributes
        MATRIX const & matrix_;
        VS vecspace_;
        magnitude_type atol_ ;
        bm_type Hdiag_ ;
        int site_ , nmin_, nmax_ ;
    };
    // -- Constructor --
    template <class MATRIX, class VS>
    davidson<MATRIX, VS>::davidson(const MATRIX& matrix, const VS& vec, const int& site, const int& nmin, const int& nmax) :
            matrix_(matrix),
            vecspace_(vec),
            site_(site),
            nmin_(nmin),
            nmax_(nmax)
    {
        vector_type tmp = new_vector(vecspace_) ;
        Hdiag_ = contraction::diagonal_hamiltonian(matrix_.left, matrix_.right, matrix_.mpo, tmp);
    } ;
    // -- Method used to compute eigenpairs --
    template <class MATRIX, class VS>
    template <class GEN, class ITER>
    std::pair<typename davidson<MATRIX, VS>::magnitude_type, typename davidson<MATRIX, VS>::vector_type>
    davidson<MATRIX, VS>::calculate_eigenvalue(const GEN& gen, ITER& iter) {
        // Initialization
        vector_type t  = new_vector(vecspace_);
        vector_type u  = new_vector(vecspace_);
        vector_type uA = new_vector(vecspace_);
        vector_type r  = new_vector(vecspace_);
        vector_set V2, VA ;
        magnitude_type theta ;
        atol_ = iter.absolute_tolerance();
        std::size_t iter_dim = 0 ;
        // Generate guess
        ietl::generate(t,gen);
        ietl::project(t,vecspace_);
        // While loop
        do {
            update_vspace(V2, VA, t, iter_dim);
            iter_dim = V2.size();
            matrix_numeric M(iter_dim, iter_dim), Mevecs(iter_dim, iter_dim);
            std::vector<magnitude_type> Mevals(iter_dim);
            for (int i = 0; i < iter_dim; ++i)
                for (int j = i; j < iter_dim; ++j)
                {
                    M(i,j) = ietl::dot(V2[i], VA[j]);
                    M(j,i) = M(i,j);
                }
            // TODO ALB Use here the same diagonalization algorithm as for JD
            boost::numeric::bindings::lapack::heevd('V', M, Mevals);
            Mevecs = M ;
	        select_eigenpair(V2, VA, Mevecs, iter_dim, u, uA) ;
            magnitude_type energy = ietl::dot(u,uA)  ;
            r = uA - energy*u ;
            theta = energy ;
            // if (|r|_2 < \epsilon) stop
            ++iter;
            if (iter.finished(ietl::two_norm(r), theta))
                break;
            precondition(r, u, theta);
            // Restarting algorithm
            t = finalize_iteration(u, r, nmax_, iter_dim, V2, VA);
        } while (true);
        // accept lambda=theta and x=u
        return std::make_pair(return_final(theta), u);
    }
}
#endif
