/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2011 by Rene Villiger <rvilliger@smile.ch>,
 *                            Prakash Dayal <prakash@comp-phys.org>,
 *                            Matthias Troyer <troyer@comp-phys.org>
 *                            Bela Bauer <bauerb@phys.ethz.ch>
 *               2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef IETL_JACOBI_H
#define IETL_JACOBI_H

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/cg.h>

#include <vector>

#include <boost/function.hpp>

namespace ietl
{
    //
    // +---------------------+
    //  JACOBI-DAVIDSON CLASS
    // +---------------------+
    // This is a general class for Davidson-type eigensolver.
    // The templates arguments are MATRIX and VS, that are usually:
    // MATRIX    : a SiteProblem object (see optimize.h for additional details -
    // VS        : a VectorSpace object, including a MPSTensor and several other vectors
    //             (for excited states orthogonalization)
    //
    // Includes the following attributes, that are common to all the Davidson eigensolvers;
    // matrix_   : the matrix representation of the operator for the site where the optimization is
    //             carried out
    // vecspace_ : the vector space where the optimization is carried out
    // Other attributes, that are specific of other type of eigensolvers (such as a shift omega or a
    // state to target) are defined in the inherited classes.
    //
    // The methods are the following:
    // contructor           : standard constructor
    // calculate_eigenvalue : method to compute eigenvalue, uses virtual functions
    // update_vspace        : virtual protected function, defined when virtual class is inherited
    // apply_operator       : virtual protected function, defined when virtual class is inherited
    // precondition         : virtual protected function, for the guess iteration of each Davidson step
    //
    template <class MATRIX, class VS, class ITER>
    class jacobi_davidson
    {
    public:
        typedef typename vectorspace_traits<VS>::scalar_type                     scalar_type;
        typedef typename vectorspace_traits<VS>::vector_type                     vector_type;
        typedef typename vector_type::bm_type 		      		                 bm_type ;
        typedef typename std::vector< bm_type >                                  vector_bm ;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type        magnitude_type;
        typedef typename std::size_t                                             size_t ;
        typedef typename std::pair<int, float>                                   couple_val ;
        typedef typename std::vector<couple_val>                                 couple_vec ;
        typedef typename std::vector<double>                                     vector_double ;
        typedef typename std::vector<vector_double>                              matrix_double ;
        typedef typename std::vector<vector_type>                                vector_space;
        typedef typename std::pair<magnitude_type, vector_type >                 pair_results ;
        typedef typename std::pair<vector_type , vector_type >                   pair_vectors ;
        typedef typename std::vector< pair_results >                             vector_pairs ;
        typedef typename std::vector< std::pair < pair_vectors, vector_type > >  result_collector ;
        typedef typename std::vector< pair_vectors >                             vector_ortho_vec ;
        jacobi_davidson(const MATRIX& matrix, VS& vec, const size_t& n_min, const size_t& n_max,
                        const size_t& max_iter, const int& nsites, const int& site1, const int& site2,
                        const double& ietl_tol, const size_t & i_gmres_guess);
        virtual ~jacobi_davidson() {};
        template <class GEN>
        vector_pairs calculate_eigenvalue(const GEN& gen, ITER& iter);
    protected:
        //
        void get_eigenvalue(std::vector<double>& eigval, std::vector<class std::vector<double> >& eigvecs, fortran_int_t dim,
                            fortran_int_t i1, fortran_int_t i2) ;
        // Virtual protected methods, to be inherited by derived classes
        virtual bool check_convergence(const vector_type& u, const vector_type& uA , const vector_type& r, const magnitude_type theta ,
                                       ITER& iter, vector_type& eigvec, magnitude_type& eigval) {};
        virtual magnitude_type compute_property(const vector_space& V, const vector_space& VA, const int& i) {} ;
        virtual vector_double generate_property(const vector_space& V, const vector_space& VA, const size_t& dim,
                                                const matrix_double& eigvecs, const vector_double& eigvals) {} ;
        virtual vector_type apply_operator (const vector_type& x) {} ;
        virtual vector_type compute_error (const vector_type& u , const vector_type& uA, magnitude_type theta) {} ;
        virtual void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                            vector_type& output, vector_type& outputA, magnitude_type& theta,
                                            matrix_double& eigvecs, vector_double& eigvals) {} ;
        virtual void print_endline(void) {} ;
        virtual void print_header_table(void) {} ;
        virtual void print_newline_table(const size_t& i , const double& er, const magnitude_type& ener, const double& ov) {} ;
        virtual void solver(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t,
                            const magnitude_type& rel_tol ) {} ;
        virtual void sort_prop(couple_vec& vector_values) {} ;
        virtual void update_vecspace(vector_space &V, vector_space &VA, const int i, vector_pairs& res) {};
        virtual void update_orthospace(void) {} ;
        virtual void update_u_and_uA(const vector_type& u, const vector_type& uA) {} ;
        // Structure used for restart
        struct lt_couple {
            inline bool operator() (const couple_val& a , const couple_val& b) {
                return (a.second < b.second) ;
            }
        };
        struct gt_couple {
            inline bool operator() (const couple_val& a , const couple_val& b) {
                return (a.second > b.second) ;
            }
        };
        // Protected attributes
        double ietl_tol_, overlap_ ;
        int nsites_, site1_, site2_ ;
        FortranMatrix<scalar_type> M ;
        MATRIX const &      matrix_ ;
        vector_ortho_vec    ortho_space_left_ , ortho_space_right_ ;
        size_t              i_gmres_guess_, i_state_, max_iter_ , n_restart_min_ , n_restart_max_, n_root_found_, n_sa_ ;
        vector_bm           Hdiag_ ;
        result_collector    u_and_uA_ ;
        vector_space        v_guess_ ;
        VS                  vecspace_ ;
    private:
        // Private method, interface to the LAPACK diagonalization routine
        void restart_jd(vector_space &V, vector_space &VA, const matrix_double& eigvec, const vector_double& eigval);
        // Struct used for final sorting
        struct lt_result {
            inline bool operator() (const pair_results& a , const pair_results& b) {
                return (a.first < b.first) ;
            }
        };
    };
    // -- Constructor --
    template <class MATRIX, class VS, class ITER>
    jacobi_davidson<MATRIX, VS, ITER>::jacobi_davidson(const MATRIX& matrix, VS& vec, const size_t& n_min, const size_t& n_max,
                                                       const size_t& max_iter, const int& nsites, const int& site1, const int& site2,
                                                       const double& ietl_tol, const size_t& i_gmres_guess ) :
        matrix_(matrix),
        vecspace_(vec),
        nsites_(nsites),
        site1_(site1),
        site2_(site2),
        M(1,1),
        max_iter_(max_iter),
        n_restart_min_(n_min),
        n_restart_max_(n_max),
        n_root_found_(0),
        overlap_(0.),
        ietl_tol_(ietl_tol),
        i_gmres_guess_(i_gmres_guess),
        i_state_(0)
    {
        // Generates guess
        n_sa_ = n_root(vec) ;
        for (size_t k = 0; k < n_sa_; k++) {
            v_guess_.push_back(new_vector(vec, k));
        }
    } ;
    // -- Calculation of eigenvalue --
    template <class MATRIX, class VS, class ITER>
    template <class GEN>
    typename jacobi_davidson<MATRIX, VS, ITER>::vector_pairs
             jacobi_davidson<MATRIX, VS, ITER>::calculate_eigenvalue(const GEN& gen, ITER& iter)
    {
        // Variable declaration
        // Scalars
        magnitude_type eigval, rel_tol, theta;
        bool converged ;
        int n_iter ;
        rel_tol = ietl_tol_ ;
        // Vectors
        vector_double props(iter.max_iterations()) ;
        vector_space  V(iter.max_iterations())  ;
        vector_space  VA(iter.max_iterations()) ;
        vector_type u, uA;
        vector_double eigvals ;
        matrix_double eigvecs ;
        vector_type  eigvec ;
        vector_pairs res ;
        // Initialization
        M.resize(iter.max_iterations(), iter.max_iterations());
        print_header_table() ;
        //
        // Main loop of the algorithm
        // --------------------------
        // Loop over all the states to be orthogonalized
        for (int k = 0 ; k < n_sa_ ; k++){
            V[0] = v_guess_[k] ;
            n_iter = 0 ;
            i_state_ = k ;
            do {
                update_vecspace(V, VA, n_iter, res);
                // Update of the M matrix and compute the eigenvalues and the eigenvectors
                for (int j = 0; j < n_iter + 1; j++)
                     for (int i = 0; i < j + 1; i++)
                        M(i, j) = ietl::dot(V[i], VA[j]);
                diagonalize_and_select(V, VA, n_iter + 1, u, uA, theta, eigvecs, eigvals);
                // Check convergence
                ++iter;
                n_iter += 1;
                vector_type r = compute_error(u, uA, theta);
                converged = check_convergence(u, uA, r, theta, iter, eigvec, eigval);
                print_newline_table(n_iter, ietl::two_norm(r), eigval, overlap_);
                if (converged) {
                    print_endline();
                    n_root_found_ += 1 ;
                    eigvec /= ietl::two_norm(eigvec) ;
                    res.push_back(std::make_pair(eigval, eigvec)) ;
                    if (k != n_sa_-1) {
                        ortho_space_left_.resize(0) ;
                        ortho_space_right_.resize(0) ;
                        update_u_and_uA(eigvec, uA) ;
                        update_orthospace() ;
                    }
                    break ;
                }
                solver(u, theta, r, V[n_iter], rel_tol);
                if (n_iter == n_restart_max_) {
                    restart_jd(V, VA, eigvecs, eigvals);
                    print_endline();
                    n_iter = n_restart_min_ - 1;
                }
            } while (true);
        }
        return res ;
    }
    // Restarting routine
    template <class MATRIX, class VS, class ITER>
    void jacobi_davidson<MATRIX, VS, ITER>::restart_jd (vector_space &V,
                                                        vector_space &VA,
                                                        const matrix_double& eigvecs,
                                                        const vector_double& eigvals)
    {
        // Variable declaration
        std::vector<couple_val> vector_values ;
        vector_double p_tmp(n_restart_max_) ;
        vector_space V_tmp(n_restart_min_), VA_tmp(n_restart_min_) ;
        vector_type tmp_V ;
        size_t idx ;
        // Rotates the properties
        p_tmp = generate_property(V, VA, n_restart_max_, eigvecs, eigvals) ;
        // Build the vector
        for (int i = 0; i < n_restart_max_ ; i++)
            vector_values.push_back(std::make_pair(i,fabs(p_tmp[i])));
        sort_prop(vector_values) ;
        // Finalization
        for (int i = 0; i < n_restart_min_; i++) {
            idx = vector_values[i].first ;
            V_tmp[i]  = eigvecs[idx][0] * V[0];
            VA_tmp[i] = eigvecs[idx][0] * VA[0];
            for (int j = 1; j < n_restart_max_; ++j) {
                V_tmp[i]  += eigvecs[idx][j] * V[j];
                VA_tmp[i] += eigvecs[idx][j] * VA[j];
            }
        }
        for (int i = 0; i < n_restart_min_ ; i++){
            V[i]  = V_tmp[i] ;
            VA[i] = VA_tmp[i];
        }
    }
    // -- Interface to the LAPACK diagonalization routine --
    template <class MATRIX, class VS, class ITER>
    void jacobi_davidson<MATRIX, VS, ITER>::get_eigenvalue(std::vector<double>& eigval, std::vector<class std::vector<double> >& eigvec,
                                                           fortran_int_t dim, fortran_int_t id_min, fortran_int_t id_max)
    {
        // Definition of all the quantities needed by the LAPACK routine
        double abstol = 1.0E-7;
        char jobz  = 'V';
        char range = 'I';
        char uplo  = 'U';
        fortran_int_t n     = dim ;
        fortran_int_t lda   = dim ;
        fortran_int_t ldz   = n   ;
        fortran_int_t lwork = 8*n ;
        fortran_int_t info;
        fortran_int_t neig  = id_max-id_min+1 ;
        double vl, vu;
        double *w = new double[neig];
        double *z = new double[neig*n];
        double *work = new double[lwork];
        fortran_int_t *iwork = new fortran_int_t[5*n];
        fortran_int_t *ifail = new fortran_int_t[n];
        // Convert the matrix from general MATRIX class to a FortranMatrix object
        FortranMatrix<scalar_type> M_(dim,dim);
        for (int i=0 ; i<dim ; i++)
            for (int j=0 ; j<=i ; j++)
                M_(j, i) = M(j, i);
        LAPACK_DSYEVX(&jobz, &range, &uplo, &n, M_.data(), &lda, &vl, &vu, &id_min, &id_max, &abstol, &neig, w, z, &ldz, work, &lwork, iwork, ifail, &info);
        for (int j = 0 ; j < neig ; j++) {
            eigval[j] = w[j];
            for (int i = 0; i < n; i++)
                eigvec[j][i] = z[i + n*j];
        }
        // Free space
        delete [] w     ;
        delete [] z     ;
        delete [] work  ;
        delete [] iwork ;
        delete [] ifail ;
    };
}

#endif
