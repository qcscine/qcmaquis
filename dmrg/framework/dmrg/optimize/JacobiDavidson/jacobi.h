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
#include <dmrg/optimize/ietl_lanczos_solver.h>

#include "dmrg/optimize/POverlap/partial_overlap.h"
#include "dmrg/optimize/utils/varianceoptimizer.h"
#include "dmrg/optimize/utils/state_prop.h"

#include "dmrg/optimize/Finalizer/finalizer.h"
#include "dmrg/optimize/CorrectionEquation/correctionequation.h"
#include "dmrg/optimize/MicroOptimizer/optimizationalgorithm.h"
#include "dmrg/optimize/Orthogonalizer/orthogonalizer.h"

namespace ietl
{
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
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    class jacobi_davidson
    {
    public:
        typedef typename CorrectionEquation<MATRIX, VS>::CorrectionEquation                  CorrectionEquation ;
        typedef typename Finalizer<MATRIX, VS>::Finalizer                                    Finalizer ;
        typedef typename Orthogonalizer<VS>::Orthogonalizer                                  Orthogonalizer ;
        typedef typename OptimizationAlgorithm<MATRIX, VS, CorrectionEquation>::OptimizationAlgorithm MicroOptimizer ;
        typedef typename vectorspace_traits<VS>::scalar_type                                 scalar_type;
        typedef typename vectorspace_traits<VS>::vector_type                                 vector_type;
        typedef typename vectorspace_traits<VS>::real_type                                   real_type;
        typedef typename std::vector<real_type>                                              vector_real ;
        typedef typename std::vector<vector_real>                                            real_space ;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type                    magnitude_type;
        typedef typename std::size_t                                                         size_t ;
        typedef typename std::pair<size_t, float>                                            couple_val ;
        typedef typename std::pair<scalar_type , scalar_type>                                extreme_val ;
        typedef typename std::vector<couple_val>                                             couple_vec ;
        typedef typename std::vector<double>                                                 vector_double ;
        typedef typename std::vector<scalar_type>                                            vector_scalar ;
        typedef typename std::vector<magnitude_type>                                         vector_magnitude ;
        typedef typename std::vector<vector_double>                                          matrix_double ;
        typedef typename std::vector<vector_type>                                            vector_space ;
        typedef typename std::pair<magnitude_type, vector_type >                             pair_results ;
        typedef typename std::vector< pair_results >                                         vector_pairs ;
        typedef typename std::vector< std::vector<vector_type> >                             result_collector ;
        typedef typename VarianceOptimizer<SymmGroup>::VarianceOptimizer                     variance_optimizer ;
        typedef typename std::vector< state_prop<VS > >                                      vector_prop ;
        // -- CONSTRUCTOR --
        jacobi_davidson(MATRIX& matrix, VS& vec, CorrectionEquation& corrector, std::unique_ptr<MicroOptimizer>& micro_iterator,
                        Finalizer& finalizer, std::unique_ptr<Orthogonalizer> & ortho, const size_t& n_min, const size_t& n_max,
                        const size_t& n_block, const double& thresh_block, const int& site1, const int& site2,
                        const std::vector<std::size_t>& order, const int& sa_alg, const int& n_lanczos,
                        const bool& do_chebychev, const magnitude_type& chebyshev_shift, const bool& do_H_squared,
                        const bool& reshuffle_variance, const bool& track_variance, const bool& is_folded,
						const double& energy_thresh);
        virtual ~jacobi_davidson() {};
        vector_pairs calculate_eigenvalue(ITER& iter);
    protected:
        // -- STANDARD METHODS --
        int get_dim(const vector_space& V) ;
        void estimate_extremes(vector_type& t) ;
        pair_results generate_eigenpair(const size_t& idx) ;
        real_type compute_variance(const vector_type& eigen, const bool& print) ;
        size_t select_among_converged() ;
        vector_type filter_chebyshev(vector_type& x) ;
        vector_type generate_eigenvector(const size_t& idx) ;
        vector_type multiply_by_operator(const vector_type& input) ;
        void initialize_vecspace() ;
        void reset_vecspace() ;
        void update_orthospace(const vector_type& u) ;
        // -- VIRTUAL METHODS --
        //  Interface to the LAPACK diagonalization routine
        //  varies by the used class
        virtual void get_eigenvalue(std::vector<double>& eigval, std::vector<std::vector<double> >& eigvecs, size_t dim);

        virtual bool check_convergence(size_t const& idx, ITER& iter) = 0;
        virtual vector_double generate_property() = 0 ;
        virtual vector_type apply_operator (const vector_type& x) = 0 ;
        virtual void print_endline() = 0 ;
        virtual void print_header_table() = 0 ;
        virtual void print_newline_table(const size_t& i, const real_type& er, const scalar_type& ener,
                                         const size_t& idx, const bool& converged) = 0 ;
        virtual void set_interval(const std::size_t& dim) = 0 ;
        virtual void solver(const vector_type& r, vector_space& t) = 0 ;
        virtual void sort_prop(couple_vec& vector_values) = 0 ;
        virtual void update_finalizer() = 0 ;
        virtual void update_parameters() = 0 ;
        virtual void update_vecspace(vector_space& to_add) = 0;

        // Restarting routine
        void diagonalize(const bool& is_low) ;
        // diagonalize() wrapper, used in different classes
        virtual void do_diagonalize()
        {
            diagonalize(true);
            diagonalize(false);
        }
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
        // Diagonalization lower and upper bounds
        struct diag_bound{
            // Attributes
            bool do_lower, do_upper ;
            std::pair<size_t , size_t> lower_bounds ;
            std::pair<size_t , size_t> upper_bounds ;
            // Constructor
            diag_bound() : do_lower(false), do_upper(false),
                           lower_bounds(std::make_pair(0,0)),
                           upper_bounds(std::make_pair(0,0)) {} ;
        } ;
        // Protected attributes
        bool                            do_chebychev_, do_H_squared_, is_folded_, reshuffle_variance_, track_variance_ ;
        struct diag_bound               bounds ;
        double                          energy_thresh_, thresh_block_ ;
        int                             n_lanczos_, sa_alg_, site1_, site2_ ;
        Finalizer&                      finalizer_ ;
        std::unique_ptr<Orthogonalizer>& orthogonalizer_ ;
        FortranMatrix<magnitude_type>   M ;
        std::vector<scalar_type>        diagonal_elements_ ;
        CorrectionEquation &      corrector_ ;
        magnitude_type                  lowest_eigen_, highest_eigen_, chebyshev_shift_, energy_ref_, atol_ ;
        MATRIX&                         matrix_ ;
        std::unique_ptr<MicroOptimizer>& micro_iterator_ ;
        result_collector                u_and_uA_ ;
        std::vector<couple_val>         vector_values_ ;
        std::vector<std::size_t>        order_ ;
        size_t                          i_homing_selected_, i_state_, n_block_, n_iter_, n_restart_min_, n_restart_max_,
                                        n_root_found_, n_sa_ ;
        vector_prop                     eigen_collector_, candidates_collector_, converged_collector_, not_converged_collector_ ;
        vector_space                    v_guess_, V_, VA_ ;
        VS                              vecspace_ ;
    private:
        void select_u_and_uA(const size_t& n_block_local) ;
        void reshuffle_variance() ;
        void restart_jd(vector_space &to_add) ;
    };
    // -- Constructor --
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::jacobi_davidson
            (MATRIX& matrix, VS& vec, CorrectionEquation& corrector, std::unique_ptr<MicroOptimizer>& micro_iterator, Finalizer& finalizer,
             std::unique_ptr<Orthogonalizer> & ortho, const size_t& n_min, const size_t& n_max, const size_t& n_block, const double& thresh_block,
             const int& site1, const int& site2, const std::vector<std::size_t>& order, const int& sa_alg, const int& n_lanczos,
             const bool& do_chebychev, const magnitude_type& chebyshev_shift, const bool& do_H_squared, const bool& reshuffle_variance,
             const bool& track_variance, const bool& is_folded, const double& energy_thresh) :
        bounds(),
        chebyshev_shift_(chebyshev_shift),
        matrix_(matrix),
        vecspace_(vec),
        site1_(site1),
        site2_(site2),
        M(1,1),
        micro_iterator_(micro_iterator),
        n_block_(n_block),
        n_iter_(0),
        n_lanczos_(n_lanczos),
        n_restart_min_(n_min),
        n_restart_max_(n_max),
        n_root_found_(0),
        order_(order),
        orthogonalizer_(ortho),
        i_homing_selected_(0),
        i_state_(0),
		    is_folded_(is_folded),
        sa_alg_(sa_alg),
        corrector_(corrector),
        finalizer_(finalizer),
        thresh_block_(thresh_block),
        do_chebychev_(do_chebychev),
        do_H_squared_(do_H_squared),
        reshuffle_variance_(reshuffle_variance),
        track_variance_(track_variance),
        V_(0),
        VA_(0),
        energy_ref_(0.),
        energy_thresh_(energy_thresh),
        lowest_eigen_(0.),
        highest_eigen_(0.)
    {
        // Set various parameters
        n_sa_  = n_root(vec) ;
        // Generates guess
        v_guess_.resize(n_sa_) ;
        for (size_t k = 0; k < n_sa_; k++)
            v_guess_[k] = vec.new_vector(k) ;
        // Eigen collector object
        eigen_collector_.reserve(n_block_) ;
        candidates_collector_.reserve(n_block_) ;
        converged_collector_.reserve(n_block_) ;
    } ;
    // -- Calculation of eigenvalue --
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::vector_pairs
             jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::calculate_eigenvalue(ITER& iter)
    {
        // Scalars
        bool converged, converged_loc ;
        magnitude_type energy ;
        size_t n_block_local ;
        // Vectors
        matrix_double   eigvecs1, eigvecs2 ;
        vector_double   eigvals1, eigvals2 ;
        vector_pairs    res ;
        vector_space    v_toadd , r ;
        vector_type     eigvec, error ;
        // Initializatin
        corrector_.update_vecspace(vecspace_) ;
        atol_ = iter.absolute_tolerance();
        M.resize(iter.max_iterations()*n_block_, iter.max_iterations()*n_block_) ;
        res.resize(n_sa_) ;
        print_header_table() ;
        V_.resize(iter.max_iterations());
        VA_.resize(iter.max_iterations());
        orthogonalizer_->set_vecspace(V_) ;
        orthogonalizer_->set_addspace(VA_) ;
        orthogonalizer_->set_diagonal(diagonal_elements_) ;
        // +--------------------------+
        //  Main loop of the algorithm
        // +--------------------------+
        // Loop over all the states to be orthogonalized
        for (size_t k = 0 ; k < n_sa_ ; k++) {
            bool starting = true ;
            i_state_ = order_[k] ;
            reset_vecspace() ;
            do {
                // Initialization/update of the vector space
                if (starting) {
                    initialize_vecspace();
                    update_parameters() ;
                    if (k == 0)
                        corrector_.update_hamiltonian(matrix_) ;
                    // update preconditioner for higher roots only if we don't use
                    // the same boundary for all states
                    bool update_prec = ((sa_alg_ == -1) || (sa_alg_ == -3)) ? (k == 0) : true;
                    if (update_prec)
                        corrector_.update_n_root(i_state_) ;
                    starting = false;
                } else {
                    update_vecspace(v_toadd);
                }
                // Preliminary operators
                r.clear();
                n_iter_++;
                // Update of the M matrix and compute the eigenvalues and the eigenvectors
                for (size_t i = 0; i < V_.size(); i++)
                    M(i, n_iter_-1) = orthogonalizer_->get_hamiltonian(V_, VA_, i , V_.size()-1);
                // Diagonalization
                set_interval(n_iter_);
                do_diagonalize();
                // Check convergence
                n_block_local = std::min(candidates_collector_.size(), static_cast<size_t>(n_block_)) ;
                select_u_and_uA(n_block_local) ;
                converged = true;
                update_finalizer() ;
                not_converged_collector_.clear() ;
                for (size_t i = 0; i < eigen_collector_.size(); i++) {
                    r.push_back(finalizer_.compute_error(i_state_, i));
                    vecspace_.project(r[i]);
                    converged_loc = check_convergence(i, iter);
                    converged = converged && converged_loc;
                    real_type error = ietl::two_norm(r[i]) ;
                    print_newline_table(n_iter_, error, finalizer_.compute_energy(i_state_, i), i, converged_loc);
                    if (converged_loc) {
                        converged_collector_.emplace_back(std::move(eigen_collector_[i]));
                        vecspace_.add_within_vec(generate_eigenvector(converged_collector_.size() - 1));
                    } else {
                        not_converged_collector_.emplace_back(std::move(eigen_collector_[i]));
                    }
                }
                // Hack to add at least one vector if the optimization has reached the end
                if (n_iter_ == iter.max_iterations()) {
                    converged = true ;
                    if (converged_collector_.size() == 0)
                        converged_collector_.emplace_back(std::move(eigen_collector_[0])) ;
                }
                // -- CHECKS HOW MANY ROOTS ARE CONVERGED --
                // If converged, all the roots are converged. For SA solver, moves to the next state, otherwise just exit
                if (converged || n_iter_ == iter.max_iterations()) {
                    print_endline();
                    if (reshuffle_variance_)
                        reshuffle_variance() ;
                    size_t i_selected = select_among_converged() ;
                    res[i_state_] = generate_eigenpair(i_selected) ;
                    if (track_variance_)
                        compute_variance(res[i_state_].second, true) ;
                    n_root_found_ += 1 ;
                    if (k != n_sa_ - 1) {
                        vecspace_.clear_projector() ;
                        if (finalizer_.get_is_si())
                            update_orthospace(converged_collector_[0].uA_) ;
                        else
                            update_orthospace(converged_collector_[0].u_) ;
                        corrector_.update_vecspace(vecspace_);
                        n_iter_ = 0 ;
                    }
                    converged_collector_.clear() ;
                    iter.reset();
                    break;
                // If only some roots are converged, restart using the approximation of the others as a starting guess
                } else if (not_converged_collector_.size() != eigen_collector_.size()) {
                    for (std::size_t idx = 0; idx < not_converged_collector_.size(); idx++)
                        v_toadd.push_back(not_converged_collector_[idx].u_);
                    starting = true ;
                    reset_vecspace() ;
                // If no roots are converged, updates corrector object and generates a new vector
                } else {
                    for (size_t i = 0; i < not_converged_collector_.size(); i++) {
                        r[i] = finalizer_.compute_error(i_state_, i);
                        vecspace_.project(r[i]);
                        corrector_.update_u(not_converged_collector_[i].u_);
                        corrector_.update_Au(not_converged_collector_[i].uA_);
                        corrector_.update_rayleigh();
                        corrector_.update_error(r[i]);
                        if (do_chebychev_ && V_.size() >= 2) {
                            v_toadd.emplace_back(filter_chebyshev(not_converged_collector_[i].u_)) ;
                        } else {
                            solver(r[i], v_toadd) ;
                        }
                    }
                    if (V_.size() >= n_restart_max_) {
                        restart_jd(v_toadd);
                        reset_vecspace() ;
                    }
                }
            } while (true);
        }
        return res ;
    }
    // +------------------------------------+
    //  Routine selecting the good eigenpair
    // +------------------------------------+
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::select_u_and_uA(const size_t& n_block_local)
    {
        // Initialization
        vector_double properties = generate_property() ;
        vector_values_.clear() ;
        eigen_collector_.clear() ;
        eigen_collector_.reserve(n_block_local) ;
        // Main body
        for (int i = 0; i < properties.size() ; i++) {
            vector_values_.push_back(std::make_pair(i, properties[i])) ;
            candidates_collector_[i].add_property(properties[i]) ;
        }
        sort_prop(vector_values_) ;
        i_homing_selected_ = vector_values_[0].first ;
        eigen_collector_.emplace_back(std::move(candidates_collector_[0]));
        // Finalization
        for (int i = 1; i < n_block_local; i++) {
            size_t idx2 = vector_values_[i].first ;
            magnitude_type ratio  = std::fabs(vector_values_[i].second/vector_values_[0].second) ;
            magnitude_type ratio2 = std::fabs(finalizer_.theta_converter(candidates_collector_[idx2].theta_)-energy_ref_) ;
            if ( (ratio > thresh_block_ && ratio2 < energy_thresh_) ) {
                if (eigen_collector_.size() < n_block_local)
                    eigen_collector_.emplace_back(std::move(candidates_collector_[idx2]));
            }
        }
        // At least one vector is added
        if (eigen_collector_.size() == 0)
            eigen_collector_.emplace_back(std::move(candidates_collector_[i_homing_selected_])) ;
    }
    // +---------------------------+
    //  Vector space initialization
    // +---------------------------+
    template <class Matrix, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<Matrix, VS, SymmGroup, ITER>::initialize_vecspace()
    {
        // Variables declaration
        vector_type t, tA ;
        // Orthogonalization
        // Leon: Commented the sa_alg_ == -1 orthogonalisation for now because of crashes/wrong results in some cases
        // Need to check with ALB
        // if (sa_alg_ != -1) {
            t = vecspace_.new_vector(i_state_) ;
            if (n_lanczos_ > 0)
                estimate_extremes(t) ;
            vecspace_.project(t) ;
            tA = apply_operator(t) ;
            vecspace_.project(tA) ;
            orthogonalizer_->normalize(t, tA) ;
            orthogonalizer_->update_diagonal(t, tA);
            V_.emplace_back(std::move(t)) ;
            VA_.emplace_back(std::move(tA)) ;
        // } else {
        //     //
        //     for (size_t i = 0; i < n_sa_; i++) {
        //         t = vecspace_.new_vector(i) ;
        //         vecspace_.project(t) ;
        //         tA = apply_operator(t) ;
        //         vecspace_.project(tA) ;
        //         orthogonalizer_->orthogonalize(t, tA) ;
        //         if (ietl::two_norm(tA) > 1.0E-10 && ietl::two_norm(t) > 1.0E-10 ) {
        //             orthogonalizer_->normalize(t, tA);
        //             orthogonalizer_->update_diagonal(t, tA);
        //             V_.push_back(t) ;
        //             VA_.push_back(tA) ;
        //         }
        //     }
        // }
    };
    // +---------------------+
    //  Computes the variance
    // +---------------------+
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::real_type
             jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::compute_variance(const vector_type& eigen, const bool& print)
    {
        real_type variance, jnk ;
        jnk      = get_energy(matrix_, eigen, i_state_, false) ;
        variance = get_energy(matrix_, eigen, i_state_, true) ;
        variance -= jnk*jnk ;
        if (print)
            std::cout << " Variance = " << sqrt(variance) << std::endl ;
        return sqrt(variance) ;
    }
    // +------------------+
    //  Restarting routine
    // +------------------+
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::restart_jd(vector_space& to_add)
    {
        // Check coherence in input data
        to_add.clear() ;
        assert (eigen_collector_.size() >= n_restart_min_) ;
        // Finalization
        for (int i = 0; i < n_restart_min_; i++)
            to_add.push_back(eigen_collector_[i].u_) ;
    }
    // +---------------+
    //  Diagonalization
    // +---------------+
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::diagonalize(const bool& is_low)
    {
        // -- Initialization --
        assert(V_.size() == VA_.size()) ;
        matrix_double eigvecs ;
        vector_double eigvals ;
        size_t dim = V_.size(), imin, imax ;
        // -- Upper diagonalization --
        if (is_low) {
            imin = bounds.lower_bounds.first;
            imax = bounds.lower_bounds.second;
        } else {
            imin = bounds.upper_bounds.first;
            imax = bounds.upper_bounds.second;
        }
        // -- Resizes the matrices --
        if (imin > 0) {
            eigvals.resize(dim);
            eigvecs.resize(dim);
            for (int i = 0; i < dim; i++)
                eigvecs[i].resize(dim);
            // -- Diagonalization --
            get_eigenvalue(eigvals, eigvecs, dim);
            // -- Finalization --
            couple_vec selector ;
            for (size_t i = 0; i < eigvals.size(); i++)
                selector.push_back(std::make_pair(i, eigvals[i])) ;
            //
            std::sort(selector.begin(), selector.end(), lt_couple()) ;
            candidates_collector_.clear() ;
            candidates_collector_.reserve(imax-imin+1) ;
            for (size_t i = imin-1; i < imax; i++) {
                size_t idx = selector[i].first ;
                vector_type tmp_u  = eigvecs[idx][0] * V_[0];
                vector_type tmp_uA = eigvecs[idx][0] * VA_[0];
                for (size_t j = 1; j < dim; ++j) {
                    tmp_u  += eigvecs[idx][j] * V_[j];
                    tmp_uA += eigvecs[idx][j] * VA_[j];
                }
                // tmp_uA = apply_operator(tmp_u) ;
                candidates_collector_.push_back(state_prop<VS>(tmp_u, tmp_uA, eigvals[idx])) ;
            }
        }
    }
    // +---------------------------------------------------+
    //   Reshuffle the states based on the variance metrix
    // +---------------------------------------------------+
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::reshuffle_variance()
    {
        // Initial dimensions
        std::size_t dim = converged_collector_.size() ;
        if (dim == 1)
            return ;
        vector_type jnk, jnk1, jnk2 ;
        ietl::FortranMatrix<double> H(dim, dim), H_squared(dim, dim) ;
        for (size_t idx1 = 0; idx1 < dim; idx1++) {
            for (size_t idx2 = 0; idx2 < dim; idx2++) {
                // H matrix
                if (finalizer_.get_is_si()) {
                    jnk1 = converged_collector_[idx1].uA_ ;
                    jnk2 = converged_collector_[idx2].uA_ ;
                } else {
                    jnk1 = converged_collector_[idx1].u_ ;
                    jnk2 = converged_collector_[idx2].u_ ;
                }
                jnk1 /= ietl::two_norm(jnk1) ;
                jnk2 /= ietl::two_norm(jnk2) ;
                ietl::mult(matrix_, jnk1, jnk, i_state_, false);
                H(idx1, idx2) = std::abs(ietl::dot(jnk2, jnk)) ;
                // H^2 matrix
                ietl::mult(matrix_, jnk1, jnk, i_state_, true);
                H_squared(idx1, idx2) = std::abs(ietl::dot(jnk2, jnk)) ;
            }
        }
        variance_optimizer var_opt(H, H_squared) ;
        var_opt.PerformOptimization() ;
        std::vector<double> coeff = var_opt.get_coeff() ;
        // Modifies the eigenvectors
        jnk  = converged_collector_[0].u_*coeff[0] ;
        jnk2 = converged_collector_[0].uA_*coeff[0] ;
        for (size_t idx1 = 1; idx1 < dim; idx1++) {
            jnk  += coeff[idx1] * converged_collector_[idx1].u_;
            jnk2 += coeff[idx1] * converged_collector_[idx1].uA_;
        }
        // Here we normalized just for safety
        magnitude_type theta_loc = converged_collector_[0].theta_ ;
        converged_collector_.resize(0) ;
        converged_collector_.push_back(state_prop<VS>(jnk, jnk2, theta_loc)) ;
    }
    // +---------------------------------+
    //  Get the dimension of the subspace
    // +---------------------------------+
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    int jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::get_dim(const vector_space& V)
    {
        return static_cast<int>(V.size()) ;
    }
    // +-----------------+
    //  UPDATE_ORTHOSPACE
    // +-----------------+
    // Update of the orthogonal space for the constrained optimization
    template <class Matrix, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<Matrix, VS, SymmGroup, ITER>::reset_vecspace()
    {
        V_.clear() ;
        VA_.clear() ;
        diagonal_elements_.clear() ;
    }
    // +-----------------+
    //  UPDATE_ORTHOSPACE
    // +-----------------+
    // Update of the orthogonal space for the constrained optimization
    template <class Matrix, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<Matrix, VS, SymmGroup, ITER>::update_orthospace(const vector_type &u)
    {
        for (size_t jcont = 0; jcont < n_root_found_; jcont++) {
            vector_type tmp = vecspace_.return_orthovec(u/ietl::two_norm(u), order_[n_root_found_], order_[jcont],
                                                        site1_, site2_) ;
            vecspace_.project(tmp) ;
            if (ietl::two_norm(tmp) > 1.0E-10)
                vecspace_.add_other_vec(tmp/ietl::two_norm(tmp)) ;
        }
    }
    // +----------------------+
    //  SELECT_AMONG_CONVERGED
    // +----------------------+
    // Selects the eigenvector with the largest value of the property (overlap for MaxO and energy for standard
    // Jacobi-Davidson algorithm)
    template <class Matrix, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson<Matrix, VS, SymmGroup, ITER>::size_t
             jacobi_davidson<Matrix, VS, SymmGroup, ITER>::select_among_converged()
    {
        size_t i_selected = 0 ;
        if (do_H_squared_) {
            real_type variance_ref = compute_variance(converged_collector_[i_selected].u_, false) ;
            for (size_t idx = 0; idx < converged_collector_.size(); idx++) {
                real_type variance = compute_variance(converged_collector_[idx].u_, false) ;
                if (variance < variance_ref) {
                    i_selected = idx;
                    variance_ref = variance ;
                }
            }
        } else {
            real_type property_ref = converged_collector_[i_selected].property_ ;
            for (size_t idx = 0; idx < converged_collector_.size(); idx++) {
                real_type property = converged_collector_[idx].property_ ;
                if (property > property_ref) {
                    i_selected = idx;
                    property_ref = property ;
                }
            }
        }
        return i_selected ;
    }
    // +-----------------+
    //  ESTIMATE_EXTREMES
    // +-----------------+
    template <class Matrix, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<Matrix, VS, SymmGroup, ITER>::estimate_extremes(vector_type& t)
    {
        int max_iter = n_lanczos_ ;
        if (max_iter > t.data().num_elements())
            max_iter = t.data().num_elements() - 1 ;
        vector_type jnk, jnk2 ;
        vector_space v, w ;
        v.push_back(t/ietl::two_norm(t)) ;
        vector_magnitude alpha_vec, beta_vec ;
        // +-----------------+
        //  LANCZOS ALGORITHM
        // +-----------------+
        // First iteration
        jnk = multiply_by_operator(v[0]) ;
        w.push_back(jnk) ;
        alpha_vec.push_back(std::real(ietl::dot(w[0], v[0]))) ;
        w[0] -= alpha_vec[0]*v[0] ;
        // Main loop
        std::size_t index = 1 ;
        do {
            magnitude_type norm = ietl::two_norm(w[index-1]) ;
            if (norm < 1.0E-5)
                break ;
            else
                beta_vec.push_back(norm) ;
            v.push_back( w[index-1] / beta_vec[index-1] ) ;
            jnk = multiply_by_operator(v[index]) ;
            w.push_back(jnk) ;
            alpha_vec.push_back(std::real(ietl::dot(w[index], v[index]))) ;
            w[index] -= alpha_vec[index]*v[index] + beta_vec[index-1]*v[index-1] ;
            index++ ;
        } while(index < max_iter) ;
        //
        auto actual_size = (int)v.size() ;
        // Actual diagonalization
        char jobz = 'V', range = 'A' ;
        int IL = 0, IU = 0, LDZ = actual_size, info;
        int nfound, lwork = 20*actual_size, liwork = 10*actual_size ;
        double  VL = 0., VU = 0. ;
        double tol = 1.0E-10 ;
        auto *eigval       = new magnitude_type[actual_size] ;
        auto *eigvec       = new magnitude_type[actual_size*actual_size] ;
        auto *supp_vec     = new int[2*actual_size] ;
        auto *work         = new double[lwork] ;
        auto *iwork        = new int[liwork] ;
        auto *alpha_double = new double[actual_size] ;
        auto *beta_double  = new double[actual_size-1] ;
        std::copy(alpha_vec.begin(), alpha_vec.end(), alpha_double) ;
        std::copy(beta_vec.begin(), beta_vec.end(), beta_double) ;
        LAPACK_DSTEVR(&jobz, &range, &actual_size, alpha_double, beta_double, &VL, &VU, &IL, &IU, &tol, &nfound,
                      eigval, eigvec, &LDZ, supp_vec, work, &lwork, iwork, &liwork, &info) ;
        // Updates parameters
        lowest_eigen_  = eigval[0] ;
        highest_eigen_ = eigval[actual_size-1] ;
        // Select the most similar eigenvalue
        magnitude_type reference_overlap = 0. ;
        vector_type vector_selected, v_tst ;
        for (std::size_t idx = 0; idx < v.size(); idx++) {
            v_tst = eigvec[idx*actual_size] * v[0] ;
            for (std::size_t idx2 = 1; idx2 < v.size(); idx2++)
                v_tst += eigvec[actual_size*idx + idx2] * v[idx2];
            magnitude_type overlap = std::real(ietl::dot(t / ietl::two_norm(t), v_tst / ietl::two_norm(v_tst))) ;
            if (overlap > reference_overlap) {
                vector_selected = v_tst ;
                reference_overlap = overlap ;
            }
        }
        t = vector_selected ;
        // Free memory
        delete [] eigval;
        delete [] eigvec;
        delete [] supp_vec;
        delete [] work;
        delete [] iwork;
        delete [] alpha_double;
        delete [] beta_double ;
    }
    // +--------------------+
    //  GENERATE_EIGENVECTOR
    // +--------------------+
    // Based on the orthogonalization algorithm, generates the approximation to the eigenvector
    template <class Matrix, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson<Matrix, VS, SymmGroup, ITER>::vector_type
             jacobi_davidson<Matrix, VS, SymmGroup, ITER>::generate_eigenvector(const size_t& idx)
    {
        // Compute the error vector
        vector_type eigvec ;
        if (finalizer_.get_is_si())
            eigvec = converged_collector_[idx].uA_;
        else
            eigvec = converged_collector_[idx].u_;
        eigvec /= ietl::two_norm(eigvec) ;
        return eigvec ;
    };
    // +-----------------------------------------------------------------+
    //  Generates the eigenvalue and eigenvector of the local Hamiltonian
    // +-----------------------------------------------------------------+
    template <class Matrix, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson<Matrix, VS, SymmGroup, ITER>::pair_results
             jacobi_davidson<Matrix, VS, SymmGroup, ITER>::generate_eigenpair(const size_t& idx)
    {
        // Compute the error vector
        vector_type buf , eigvec = generate_eigenvector(idx);
        magnitude_type eigval ;
        ietl::mult(this->matrix_ , eigvec , buf, i_state_, false) ;
        eigval   = std::real(ietl::dot(eigvec, buf)) ;
        return std::make_pair(eigval, eigvec) ;
    };
    // +---------------------------------------------------------------------+
    //  Filter the approximation of the eigenfunction with a Chebyshev filter
    // +---------------------------------------------------------------------+
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::vector_type
             jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::filter_chebyshev(vector_type& x)
    {
        // -- VARIABLES DECLARATION --
        vector_type x1, x2, res, x0 = x/ietl::two_norm(x) ;
        scalar_type coeff ;
        x1 = multiply_by_operator(x0) ;
        magnitude_type lower_bound = std::real(ietl::dot(x1, x0))-chebyshev_shift_;
        magnitude_type upper_bound = std::real(ietl::dot(x1, x0))+chebyshev_shift_ ;
        magnitude_type a = lowest_eigen_ ;
        magnitude_type b = highest_eigen_ ;
        magnitude_type alpha = 2./(b-a) ;
        magnitude_type beta = (b+a)/(a-b) ;
        // Check if the parameters make sense
        if (lower_bound < lowest_eigen_)
            lower_bound = lowest_eigen_ ;
        if (upper_bound > highest_eigen_)
            upper_bound = highest_eigen_ ;
        // -- ITERATIVE FILTERING --
        // First iteration
        res  = (acos(alpha*upper_bound + beta) - acos(alpha*lower_bound + beta)) * x0 / M_PI ;
        x1 *= alpha ;
        x1 += beta*x0 ;
        res += 2. * (sin( acos(alpha*upper_bound + beta) ) - sin( acos(alpha*lower_bound + beta) ) ) * x1 / M_PI ;
        // -- Main loop --
        for (std::size_t idx = 2; idx < 20; idx++) {
            // Variable setting
            coeff = 2. * (sin( idx * acos(alpha*upper_bound + beta) )
                        - sin( idx * acos(alpha*lower_bound + beta) ) ) / (M_PI*idx) ;
            // T_{n+1}(alpha*x+beta) = 2*x*T_n(alpha*x+b)
            ietl::mult(this->matrix_, x1, x2, i_state_, false) ;
            x2 *= alpha ;
            x2 += beta*x1 ;
            x2 *= 2. ;
            // T_{n+1}(alpha*x+beta) = 2*x*T_n(alpha*x+b) - T_{n-1}(alpha*x+b)
            x2 -= x0 ;
            res += coeff*x2 ;
            // Final update of the vectors
            x0    = x1 ;
            x1    = x2 ;
        }
        return res/ietl::two_norm(res) ;
    };
    // +--------------------+
    //  MULTIPLY_BY_OPERATOR
    // +--------------------+
    // Simple wrapper for the multiplication by a matrix
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::vector_type
             jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::multiply_by_operator(const vector_type& input)
    {
        vector_type output ;
        ietl::mult(this->matrix_, input, output, i_state_, false) ;
        return output ;
    }
    // +-----------------------------------------------+
    //  Interface to the LAPACK diagonalization routine
    // +-----------------------------------------------+
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson<MATRIX, VS, SymmGroup, ITER>::get_eigenvalue(std::vector<double>& eigval,
                                                                      std::vector<class std::vector<double> >& eigvecs,
                                                                      std::size_t dim)
    {
        // Definition of all the quantities needed by the LAPACK routine
        char jobvl  = 'N';
        char jobvr  = 'V';
        auto n = static_cast<fortran_int_t>(dim) ;
        fortran_int_t lda   = n ;
        fortran_int_t ldz1  = n ;
        fortran_int_t ldz2  = n ;
        fortran_int_t lwork = 8*n ;
        fortran_int_t info;
        auto *iwork = new fortran_int_t[5*n];
        auto *ifail = new fortran_int_t[n];
        auto *w1 = new magnitude_type[n];
        auto *w2 = new magnitude_type[n];
        auto *z1 = new magnitude_type[n*n];
        auto *z2 = new magnitude_type[n*n];
        auto *work = new magnitude_type[lwork];
        // Convert the matrix from general MATRIX class to a FortranMatrix object
        FortranMatrix<magnitude_type> M_(dim,dim);
        for (std::size_t i = 0 ; i<dim ; i++)
            for (std::size_t j = 0 ; j<=i ; j++)
                M_(j, i) = M(j, i);
        LAPACK_DGEEV(&jobvl, &jobvr, &n, M_.data(), &lda, w1, w2, z1, &ldz1, z2, &ldz2, work, &lwork, &info);
        for (int j = 0 ; j < dim ; j++) {
            eigval[j] = w1[j];
            for (int i = 0; i < n; i++)
                eigvecs[j][i] = z2[i+n*j];
        }
        // Free space
        delete [] w1, w2 ;
        delete [] z1, z2 ;
        delete [] work  ;
        delete [] iwork ;
        delete [] ifail ;
    }

}

// -- Derived classes --

#include "dmrg/optimize/JacobiDavidson/jacobi_standard.h"
#include "dmrg/optimize/JacobiDavidson/jacobi_modified.h"
#include "dmrg/optimize/JacobiDavidson/jacobi_standard_mo.h"
#include "dmrg/optimize/JacobiDavidson/jacobi_modified_mo.h"

#endif
