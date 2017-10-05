/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef IETL_DAVIDSON_STANDARD_MO_H
#define IETL_DAVIDSON_STANDARD_MO_H

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/cg.h>
#include <vector>

#include "dmrg/optimize/davidson.h"
#include "dmrg/optimize/utils/orthogonalizers.hpp"
#include "dmrg/optimize/partial_overlap.h"

namespace ietl {
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    class davidson_standard_mo : public davidson<MATRIX, VS> {
    public:
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap ;
        typedef typename std::vector<partial_overlap> pov_vec ;
        typedef davidson<MATRIX, VS> base;
        typedef typename base::bm_type               bm_type;
        typedef typename base::matrix_numeric        matrix_numeric;
        typedef typename base::magnitude_type        magnitude_type;
        typedef typename base::size_t                size_t;
        typedef typename base::vector_numeric        vector_numeric ;
        typedef typename base::vector_ortho_vec      vector_ortho_vec ;
        typedef typename base::vector_set            vector_set;
        typedef typename base::vector_type           vector_type;
        // Inherited attributes
        using base::Hdiag_;
	    using base::i_state_ ;
        using base::matrix_;
        using base::nsites_ ;
        using base::n_root_found_ ;
        using base::n_sa_ ;
        using base::order_ ;
        using base::ortho_space_ ;
        using base::overlap_ ;
        using base::site1_ ;
        using base::site2_ ;
        using base::printer_ ;
        using base::u_and_uA_ ;
        using base::vecspace_ ;
        using base::v_guess_ ;
        // New constructors
        davidson_standard_mo(const MATRIX &matrix, VS &vec, const pov_vec poverlap, const int& nmin, const int& nmax,
                             const int& nsites, const int& site1, const int& site2, const std::vector<int>& order,
                             const int& root_homing_type)
                : base::davidson(matrix, vec, nmin, nmax, nsites, site1, site2, order),
                  pov_(poverlap), root_homing_type_(root_homing_type)
        {};
        ~davidson_standard_mo() {};
    private:
        // Private methods
        double compute_overlap(const vector_type& vec_test) ;
        magnitude_type return_final(const magnitude_type &x) { return x; };
        vector_type apply_operator(const vector_type &x);
        vector_type finalize_iteration(const vector_type& u, const vector_type& r, const size_t& n_restart,
                                       size_t& iter_dim, vector_set& V2, vector_set& VA);
        vector_type compute_error (const vector_type& u , const vector_type& uA, magnitude_type theta) ;
        void precondition(vector_type &r, const vector_type &V, const vector_type &VA, const magnitude_type &theta,
                          const size_t& idx);
	    void select_eigenpair(const vector_set& V, const vector_set& VA, const matrix_numeric& eigvecs,
                              const vector_numeric& eigvals, const size_t& i, vector_type& u, vector_type& uA,
                              magnitude_type& theta);
        void update_vspace(vector_set &V, vector_set &VA, vector_set &t);
        void update_u_and_uA(const vector_type& u, const vector_type& uA) ;
        void update_orthospace(void) ;
        // Printing-related methods
        void print_header_table(void) ;
        void print_endline(void) ;
        void print_newline_table(const size_t& iter, const size_t& size, const magnitude_type& error, const magnitude_type& energy) ;
        // Private attribute
        int root_homing_type_ ;
        pov_vec pov_ ;
    };
    // Definition of the virtual function update_vspace
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::update_vspace(vector_set &V, vector_set &VA, vector_set &t)
    {
        size_t n_lin ;
        for (typename vector_ortho_vec::iterator it = ortho_space_.begin(); it != ortho_space_.end(); it++)
            t[0] -= ietl::dot((*it)[0], t[0]) * (*it)[0]  ;
        n_lin = gram_schmidt_orthogonalizer<vector_type, magnitude_type>(V, t) ;
        assert (V.size()-VA.size() == n_lin) ;
        for (typename vector_set::iterator it = V.begin()+VA.size() ; it != V.end() ; it++)
            VA.push_back(apply_operator(*it));
    };
    // Definition of the virtual function apply_operator
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    typename davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::vector_type
             davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::apply_operator(const vector_type &x) {
        vector_type tmp;
        ietl::mult(matrix_, x, tmp, i_state_);
        return tmp;
    };
    // Compute the error vector
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    typename davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::vector_type
             davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::compute_error(const vector_type &u,
                                                                                     const vector_type &uA,
                                                                                     magnitude_type theta)
    {
        vector_type r = uA ;
        r -= theta*u;
        // Deflates the error vector
        for (typename vector_ortho_vec::iterator it = ortho_space_.begin(); it != ortho_space_.end(); it++)
            if (ietl::dot((*it)[0], (*it)[0]) > 1.0E-15)
                r -= ietl::dot((*it)[0],r) * (*it)[0] ;
        return r ;
    }
    // Routine doing deflation
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::update_orthospace(void)
    {
        for (size_t jcont = 0; jcont < n_root_found_; jcont++) {
            vector_type tmp = vecspace_.return_orthovec(u_and_uA_[jcont][0], order_[n_root_found_], order_[jcont], site1_, site2_) ;
            for (size_t j = 0 ; j < ortho_space_.size() ; j++)
                tmp -= ietl::dot(ortho_space_[j][0], tmp) * ortho_space_[j][0] ;
            if (ietl::two_norm(tmp) > 1.0E-20) {
                tmp /= ietl::two_norm(tmp);
                std::vector< vector_type > junk ;
                junk.push_back(tmp) ;
                ortho_space_.push_back(junk);
            }
        }
    }
    // Update the vector with the quantity to orthogonalize
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::update_u_and_uA(const vector_type &u, const vector_type &uA)
    {
        vector_type tmp = u / ietl::two_norm(u) ;
        std::vector< vector_type > junk ;
        junk.push_back(tmp) ;
        u_and_uA_.push_back(junk) ;
    }
    // Definition of the virtual function precondition
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::precondition(vector_type &r, const vector_type &V, const vector_type& VA,
                                                                                const magnitude_type &theta, const size_t& idx)
    {
        magnitude_type denom, x2, x1 = ietl::dot(V, r) ;
        vector_type Vcpy = r - V * x1;
        bm_type &data = Vcpy.data();
        assert(shape_equal(data, Hdiag_[idx]));
        for (size_t b = 0; b < data.n_blocks(); ++b) {
            for (size_t i = 0; i < num_rows(data[b]); ++i) {
                for (size_t j = 0; j < num_cols(data[b]); ++j) {
                    denom = Hdiag_[idx][b](i, j) - theta;
                    if (std::abs(denom))
                        data[b](i, j) /= denom;
                }
            }
        }
        x2 = ietl::dot(V, Vcpy);
        r = Vcpy - x2 * V;
    };
    // Virtual function finalize_iteration
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    typename davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::vector_type
             davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::finalize_iteration
            (const vector_type &u, const vector_type &r, const size_t &n_restart,
            size_t &iter_dim, vector_set &V2, vector_set &VA)
    {
        vector_type result ;
        result = r ;
        return result  ;
    }
    // Routine to select the proper eigenpair
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::select_eigenpair(const vector_set& V, 
                                                                                    const vector_set& VA, 
                                                                                    const matrix_numeric& Mevecs,
										                                            const vector_numeric& Mevals,
                                                                                    const size_t& dim, 
                                                                                    vector_type& u,
                                                                                    vector_type& uA,
										                                            magnitude_type& theta)
    {
        int idx ;
        double scr ;
        std::vector<double> overlaps(dim) ;
        vector_type u_local ;
        // Compute the various overlaps
        idx = 0 ;
        for (size_t i = 0; i < dim; ++i) {
            // Conversion to the original basis
            u_local = Mevecs(0, i) * V[0];
            for (size_t j = 1; j < dim; ++j)
                u_local += Mevecs(j, i) * V[j];
            overlaps[i] = compute_overlap(u_local);
        }
        for (size_t i = idx; i < dim; ++i) {
            if (overlaps[i] > overlaps[idx])
                idx = i;
        }
        overlap_ = overlaps[idx] ;
	    // Finalization
        u  = V[0]*Mevecs(0,idx) ;
        uA = VA[0]*Mevecs(0,idx);
        for (int i = 1; i < dim; ++i) {
            u  += Mevecs(i, idx) * V[i];
            uA += Mevecs(i, idx) * VA[i];
        }
        uA /= ietl::two_norm(u) ;
        u  /= ietl::two_norm(u) ;
	    theta = Mevals[idx] ;
        return ;
    }
    // Routine to compute the overlaps
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    double davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::compute_overlap(const vector_type &vec_test)
    {
        double ret1, ret2, ret ;
        // Calculates the two overlaps
        if (root_homing_type_ == 1 || root_homing_type_ == 3) {
            if (nsites_ == 1)
                ret1 = pov_[i_state_].overlap(vec_test/ietl::two_norm(vec_test), site1_);
            else
                ret1 = pov_[i_state_].overlap(vec_test/ietl::two_norm(vec_test), site1_, site2_);
        }
        if (root_homing_type_ == 2 || root_homing_type_ == 3)
            ret2 = ietl::dot(vec_test, v_guess_[i_state_]) / (ietl::two_norm(vec_test)*ietl::two_norm(v_guess_[i_state_]));
        // Finalizes the calculation
        if (root_homing_type_ == 1)
            ret = fabs(ret1) ;
        else if (root_homing_type_ == 2)
            ret = fabs(ret2) ;
        else
            ret = (fabs(ret1) + fabs(ret2)) * 0.5 ;
        return ret ;
    }
    // Routine to print the header of the table
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::print_header_table(void) {
        printer_.print_header_table_overlap() ;
    } ;
    //
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::print_endline(void) {
        printer_.print_endline_overlap() ;
    } ;
    //
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::print_newline_table(const size_t& iter,
                                                                                       const size_t& size,
                                                                                       const magnitude_type& error,
                                                                                       const magnitude_type& energy)
    {
        printer_.print_newline_table_overlap(iter, size, error, energy, overlap_) ;
    } ;
    //
}

#endif
