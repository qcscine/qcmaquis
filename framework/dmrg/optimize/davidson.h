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

#include <vector>

// +--------------+
//  DAVIDSON CLASS
// +--------------+
// This is a general class for Davidson-type eigensolver.
// The templates arguments are MATRIX and VS, that are usually:
// MATRIX    : a SiteProblem object (see optimize.h for additional details)
// VS        : a VectorSpace object, including a MPSTensor (or TwoSiteTensor for two site optimization)
//             and several other vectors (for excited states orthogonalization)
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
        typedef typename std::vector< bm_type >                           vector_bm ;
        typedef typename std::vector< bool >                              vector_bool ;
        typedef typename std::vector<vector_type> 	     		          vector_set ;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        typedef typename std::vector< magnitude_type >                    vector_magnitude ;
        typedef typename alps::numeric::matrix<scalar_type> 		      matrix_numeric ;
        typedef typename std::size_t 					                  size_t ;
        typedef typename std::pair<magnitude_type, vector_type >          pair_results ;
        typedef typename std::vector< pair_results >                      vector_pairs ;
        // Constructor and destructor
        davidson(const MATRIX& matrix, const VS& vec, const int& nmin, const int& nmax,
                 const int& nsites, const int& site1, const int& site2);
        virtual ~davidson() {};
        // Public method to compute eigenvalue
        template <class ITER>
        vector_pairs calculate_eigenvalue(ITER& iter);
    protected:
        // Virtual Methods, defined in the derived classes
        virtual magnitude_type return_final(const magnitude_type& eigval) {} ;
        virtual vector_type finalize_iteration(const vector_type& u, const vector_type& r, const size_t& n_restart,
                                              size_t& iter_dim, vector_set& v2, vector_set& VA) {} ;
        virtual vector_type apply_operator(const vector_type& x) {} ;
        virtual void precondition(vector_type& r, const vector_type& V, const vector_type& VA, const magnitude_type& theta,
                                  const size_t& idx) {} ;
	    virtual size_t select_eigenpair(const vector_set& V, const vector_set& VA, const matrix_numeric& eigvecs,
                                        const size_t& i, vector_set& u, vector_set& uA) {} ;
        virtual void update_vspace(vector_set& V, vector_set& VA, vector_set& t) {} ;
        // Attributes
        vector_bm Hdiag_ ;
        int site1_ , site2_, nmin_, nmax_, nsites_ ;
        size_t n_sa_  ;
        magnitude_type atol_ ;
        MATRIX const & matrix_;
        vector_bool convergence_check_ ;
        vector_set v_guess_ ;
        VS vecspace_;
    };
    // -- Constructor --
    template <class MATRIX, class VS>
    davidson<MATRIX, VS>::davidson(const MATRIX& matrix, const VS& vec, const int& nmin, const int& nmax,
                                   const int& nsites, const int& site1, const int& site2) :
            matrix_(matrix),
            vecspace_(vec),
            nmin_(nmin),
            nmax_(nmax),
            site1_(site1),
            site2_(site2),
            nsites_(nsites)
    {
        n_sa_ = n_root(vec) ;
        vector_type tst = new_vector(vec) ;
        for (size_t k = 0 ; k < n_sa_; k++) {
            v_guess_.push_back(new_vector_sa(vec, k));
            convergence_check_.push_back(false) ;
            Hdiag_.push_back(contraction::diagonal_hamiltonian(matrix_.left, matrix_.right, matrix_.mpo, v_guess_[k])) ;
        }
    } ;
    // -- Method used to compute eigenpairs --
    template <class MATRIX, class VS>
    template <class ITER>
    typename davidson<MATRIX, VS>::vector_pairs davidson<MATRIX, VS>::calculate_eigenvalue(ITER& iter)
    {
        // Initialization
        size_t iter_dim = 0 ;
        vector_pairs res ;
        vector_set t = v_guess_ ;
        vector_set u(n_sa_), uA(n_sa_), r(n_sa_) ;
        vector_set V, VA ;
        vector_magnitude theta(n_sa_) ;
        atol_ = iter.absolute_tolerance();
        // Generate guess
        for (size_t i = 0 ; i < n_sa_ ; i++)
            ietl::project(t[i],vecspace_);
        // While loop
        do {
            // Generate guess
            update_vspace(V, VA, t);
            iter_dim = V.size();
            matrix_numeric M(iter_dim, iter_dim), Mevecs(iter_dim, iter_dim);
            std::vector<magnitude_type> Mevals(iter_dim);
            for (int i = 0; i < iter_dim; ++i)
                for (int j = i; j < iter_dim; ++j)
                {
                    M(i,j) = ietl::dot(V[i], VA[j]);
                    M(j,i) = M(i,j);
                }
            // TODO ALB Use here the same diagonalization algorithm as for JD
            boost::numeric::bindings::lapack::heevd('V', M, Mevals);
            Mevecs = M ;
	        size_t n_chosen = select_eigenpair(V, VA, Mevecs, iter_dim, u, uA) ;
            ++iter;
            std::cout << iter_dim << std::endl ;
            for (size_t i = 0 ; i < n_chosen ; i++){
                theta[i] = ietl::dot(u[i],uA[i])/ietl::dot(u[i],u[i])  ;
                r[i]     = uA[i] - theta[i]*u[i] ;
                std::cout << theta[i] << std::endl ;
                std::cout << ietl::two_norm(r[i]) << std::endl ;
                // if (|r|_2 < \epsilon) stop
                if (iter.finished(ietl::two_norm(r[i]), theta[i]))
                    convergence_check_[i] = true;
                else
                    precondition(r[i], u[i], uA[i], theta[i], i);
            }
            // Restarting algorithm
            size_t n_sa_local = 0 ;
            for (size_t k = 0 ; k < n_chosen ; k++) {
                if (not convergence_check_[k]) {
                    t[n_sa_local] = finalize_iteration(u[k], r[k], nmax_, iter_dim, V, VA);
                    n_sa_local += 1;
                }
            }
            if (n_sa_local == 0) {
                for (size_t k = 0 ; k < n_sa_ ; k++)
                    res.push_back(std::make_pair(return_final(theta[k]), u[k])) ;
                break;
            }
            std::cout << "Cycle" << std::endl ;
        } while(true);
        // accept lambda=theta and x=u
        std::cout << "Going out" << std::endl ;
        return res ;
    }
}
#endif
