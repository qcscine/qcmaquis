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

#ifndef IETL_JACOBI_STANDARD_MO_H
#define IETL_JACOBI_STANDARD_MO_H

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/cg.h>
#include <vector>
#include "jacobi_standard.h"

// +-------------------------------------------------+
//  JACOBI-DAVIDSON EIGENSOLVER WITH OVERLAP TRACKING
// +-------------------------------------------------+
namespace ietl
{
    template <class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    class jacobi_davidson_standard_mo : public jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>
    {
    public:
        typedef jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER> base;
        typedef jacobi_davidson<MATRIX, VS, SymmGroup, ITER> grandparent;
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap;
        typedef typename std::vector<partial_overlap> pov_vec_type ;
        typedef typename base::CorrectionEquation     CorrectionEquation ;
        typedef typename base::couple_vec             couple_vec ;
        typedef typename base::Finalizer              Finalizer ;
        typedef typename base::gt_couple              gt_couple ;
        typedef typename base::lt_couple              lt_couple ;
        typedef typename base::magnitude_type         magnitude_type;
        typedef typename base::MicroOptimizer         MicroOptimizer;
        typedef typename base::Orthogonalizer         Orthogonalizer;
        typedef typename base::real_type              real_type ;
        typedef typename base::scalar_type            scalar_type ;
        typedef typename base::vector_double          vector_double ;
        typedef typename base::vector_space           vector_space ;
        typedef typename base::vector_type            vector_type ;
        //
        using base::bounds ;
        using base::candidates_collector_ ;
        using base::compute_variance ;
        using base::eigen_collector_ ;
        using base::i_homing_selected_ ;
        using base::i_state_ ;
        using base::vecspace_ ;
        using base::v_guess_ ;
        using base::V_ ;
        using grandparent::get_eigenvalue;
        // TODO: make sure this is correct and here the grandparent's get_eigenvalue() is really called, and make sure
        // it's really needed
        using grandparent::do_diagonalize;
        // same for do_diagonalize
        jacobi_davidson_standard_mo(MATRIX& matrix, VS& vec, CorrectionEquation& corrector, std::shared_ptr<MicroOptimizer>& micro_iterator,
                                    Finalizer& finalizer, std::shared_ptr<Orthogonalizer> & ortho, const pov_vec_type& pov, const int& n,
                                    const int& side_tofollow, const size_t& nmin, const size_t& nmax, const size_t& n_block,
                                    const double& block_thresh, const int& site1, const int& site2, const std::vector<std::size_t>& order,
                                    const int& sa_alg, const int& n_lanczos, const bool& do_chebychev, const magnitude_type& chebyshev_shift,
                                    const bool& do_H_squared, const bool& reshuffle_variance, const bool& track_variance,
                                    const bool& is_folded, const double& energy_thresh, const int& root_homing_type,
									const float& homing_ratio)
                : base::jacobi_davidson_standard(matrix, vec, corrector, micro_iterator, finalizer, ortho, nmin, nmax,
                                                 n_block, block_thresh, site1, site2, order, sa_alg, n_lanczos, do_chebychev,
                                                 chebyshev_shift, do_H_squared, reshuffle_variance, track_variance, is_folded,
                                                 energy_thresh) ,
                  pov_(pov) , n_maxov_(n), root_homing_type_(root_homing_type), side_tofollow_(side_tofollow),
                  homing_ratio_(homing_ratio)
        {
            // Set the bound object
            if (side_tofollow_ == 0) {
                bounds.do_lower = true ;
                bounds.do_upper = true ;
            } else if (side_tofollow_ == 1 ) {
                bounds.do_lower = true ;
            } else if (side_tofollow_ == -1 ) {
                bounds.do_upper = true ;
            } else {
                throw std::runtime_error("Unrecognized value for extreme to be diagonalized") ;
            }
            // Check the ratio
            if (homing_ratio < 0. || homing_ratio > 1.)
                throw std::runtime_error("ua_homing_ratio must be > 0 and < 1") ;
        } ;
        ~jacobi_davidson_standard_mo() {} ;
    private:
        double compute_overlap(const vector_type& u, const vector_type& uA) ;
        vector_double generate_property() ;
        void print_endline() ;
        void print_header_table() ;
        void print_newline_table(const size_t& i, const real_type& error, const scalar_type& en,
                                 const size_t& idx, const bool& converged) ;
        void set_interval(const std::size_t& dim) ;
        void sort_prop(couple_vec& vector_values) ;
        // Private attributes
        int root_homing_type_ , side_tofollow_ ;
        double homing_ratio_ ;
        pov_vec_type pov_  ;
        int n_maxov_  ;
    };
    // -- GENERATE_PROPERTY --
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    typename jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::vector_double
             jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::generate_property()
    {
        // Variable declaration
        std::size_t dim = candidates_collector_.size() ;
        vector_type tmp_V, tmp_VA ;
        vector_double p_tmp(dim) ;
        // Rotates the properties
        for (int i = 0; i < dim ; i++) {
            tmp_V    = candidates_collector_[i].u_ ;
            tmp_VA   = candidates_collector_[i].uA_ ;
            p_tmp[i] = compute_overlap(tmp_V, tmp_VA) ;
        }
        return p_tmp ;
    }
    // -- COMPUTE_OVERLAP --
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    double jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::compute_overlap(const vector_type& u,
                                                                                                  const vector_type& uA)
    {
        double ret;
        vector_type tmp_u = u/ietl::two_norm(u) ;
        vector_type tmp_uA = u/ietl::two_norm(uA) ;
        // Calculates the two overlaps
        if (root_homing_type_ == 4) {
            ret = compute_variance(tmp_u, false);
        } else {
            double ret1u, ret1uA, ret2u, ret2uA;
            if (root_homing_type_ == 1 || root_homing_type_ == 3) {
                ret1u  = fabs(pov_[i_state_].overlap(tmp_u)) ;
                ret1uA = fabs(pov_[i_state_].overlap(tmp_uA)) ;
            }
            if (root_homing_type_ == 2 || root_homing_type_ == 3) {
                ret2u  = fabs(ietl::dot(u, v_guess_[i_state_])) / (ietl::two_norm(u) * ietl::two_norm(v_guess_[i_state_]));
                ret2uA = fabs(ietl::dot(uA, v_guess_[i_state_])) / (ietl::two_norm(uA) * ietl::two_norm(v_guess_[i_state_]));
            }
            // Finalizes the calculation
            if (root_homing_type_ == 1)
                ret = (1. - homing_ratio_) * fabs(ret1u) + homing_ratio_ * fabs(ret1uA);
            else if (root_homing_type_ == 2)
                ret = (1. - homing_ratio_) * fabs(ret2u) + homing_ratio_ * fabs(ret2uA);
            else
                ret = ((1. - homing_ratio_) * fabs(ret1u) + homing_ratio_ * fabs(ret1uA)
                       + (1. - homing_ratio_) * fabs(ret2u) + homing_ratio_ * fabs(ret2uA)) * 0.5;
        }
        return ret;
    }
    // -- SET_INTERVAL --
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::set_interval(const std::size_t& dim)
    {
        // Both sides
        if (bounds.do_lower && bounds.do_upper) {
            if (dim <= 2*n_maxov_) {
                bounds.lower_bounds = std::make_pair(1, dim) ;
                bounds.upper_bounds = std::make_pair(0, 0) ;
            } else {
                bounds.lower_bounds = std::make_pair(1, n_maxov_);
                bounds.upper_bounds = std::make_pair(dim-n_maxov_+1, dim) ;
            }
        // Only lower side
        } else if (bounds.do_lower) {
            if (dim < n_maxov_)
                bounds.lower_bounds = std::make_pair(1, dim) ;
            else
                bounds.lower_bounds = std::make_pair(1, n_maxov_) ;
        // Only upper side
        } else if (bounds.do_upper) {
            if (dim < n_maxov_)
                bounds.upper_bounds = std::make_pair(1, dim) ;
            else
                bounds.upper_bounds = std::make_pair(dim-n_maxov_+1, dim) ;
        }
    }
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::sort_prop(couple_vec& vector_values)
    {
        std::sort(vector_values.begin(),vector_values.end(),gt_couple()) ;
    } ;
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_header_table()
    {
        print_endline() ;
        std::cout << " Iteration | Sub. Dim. |    Error    |    Energy    |    Property   | Hom. Indx." << std::endl ;
        print_endline() ;
    } ;
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_endline()
    {
        std::cout << "-----------+-----------+-------------+--------------+---------------+-------------" << std::endl ;
    } ;
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_newline_table
            (const size_t& i, const real_type& error, const scalar_type& en, const size_t& idx, const bool& converged)
    {
        char buf[100] ;
        int a = i, n;
        n = sprintf(buf, "%5d      |%5d      | %1.4E  |  %6.5f  |    %4.4f     | %5d", a, this->get_dim(V_), error, en,
                    eigen_collector_[idx].property_, static_cast<int>(i_homing_selected_));
        std::cout << buf ;
        if (converged)
            std::cout << " CONVERGED" ;
        std::cout << std::endl ;
    }

}

#endif
