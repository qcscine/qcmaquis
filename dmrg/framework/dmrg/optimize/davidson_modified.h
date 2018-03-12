/*****************************************************************************
 *
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

#ifndef IETL_DAVIDSON_MODIFIED_H
#define IETL_DAVIDSON_MODIFIED_H

#include <ietl/cg.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/traits.h>
#include <vector>

#include "dmrg/optimize/davidson.h"
#include "dmrg/optimize/utils/orthogonalizers.hpp"

namespace ietl {
    template<class MATRIX, class VS>
    class davidson_modified : public davidson<MATRIX, VS> {
    public:
        typedef davidson<MATRIX, VS> base;
        typedef typename base::bm_type               bm_type ;
        typedef typename base::magnitude_type        magnitude_type ;
    	typedef typename base::matrix_numeric        matrix_numeric ;
        typedef typename base::size_t                size_t ;
        typedef typename base::vector_numeric        vector_numeric ;
        typedef typename base::vector_ortho_vec      vector_ortho_vec ;
        typedef typename base::vector_set            vector_set ;
        typedef typename base::vector_type           vector_type ;
        // Inherited attributes
        using base::Hdiag_ ;
        using base::i_state_ ;
        using base::matrix_ ;
        using base::n_root_found_ ;
        using base::n_sa_ ;
        using base::nsites_ ;
        using base::order_ ;
        using base::ortho_space_ ;
        using base::printer_ ;
        using base::site1_ ;
        using base::site2_ ;
        using base::u_and_uA_ ;
        using base::vecspace_ ;
        // New constructors
        davidson_modified(const MATRIX &matrix, VS &vec, const magnitude_type& omega, const int& nmin,
                          const int& nmax, const int& nsites, const int& site1, const int& site2,
                          const std::vector<int>& order)
                : base::davidson(matrix, vec, nmin, nmax, nsites, site1, site2, order) , omega_(omega) {};
        ~davidson_modified() {};
    private:
        // Private methods
        magnitude_type return_final(const magnitude_type &x) { return omega_-1./x; };
        void select_eigenpair(const vector_set& V, const vector_set& VA, const matrix_numeric& eigvecs,
                              const vector_numeric& eigvals, const size_t& i, vector_type& u, vector_type& uA,
                              magnitude_type& theta);
        vector_type apply_operator(const vector_type &x);
        vector_type compute_error (const vector_type& u , const vector_type& uA, magnitude_type theta) ;
        vector_type finalize_iteration(const vector_type& u, const vector_type& r, const size_t& n_restart,
                                       size_t& iter_dim, vector_set& V2, vector_set& VA);
        void precondition(vector_type &r, const vector_type &V, const vector_type &VA, const magnitude_type &theta,
                          const size_t& idx);
        void update_vspace(vector_set &V, vector_set &VA, vector_set &t);
        void update_u_and_uA(const vector_type& u, const vector_type& uA) ;
        void update_orthospace(void) ;
        // Printing-related methods
        void print_header_table(void) ;
        void print_endline(void) ;
        void print_newline_table(const size_t& iter, const size_t& size, const magnitude_type& error, const magnitude_type& energy) ;
        // Additional attributes
        magnitude_type omega_ ;
        vector_set V_additional_ ;
    };
    // Construction of the virtual function apply operator
    template <class MATRIX, class VS>
    typename davidson_modified<MATRIX, VS>::vector_type davidson_modified<MATRIX, VS>::apply_operator(const vector_type& x)
    {
        vector_type tmp ;
        ietl::mult(matrix_, x, tmp, i_state_);
        tmp *= -1. ;
        tmp += omega_*x ;
        return tmp ;
    };
    // Compute the error vector
    template <class MATRIX, class VS>
    typename davidson_modified<MATRIX, VS>::vector_type davidson_modified<MATRIX,VS>::compute_error(const vector_type &u,
                                                                                                    const vector_type &uA,
                                                                                                    magnitude_type theta)
    {
        vector_type r ;
        r = uA - u/theta ;
        r /= ietl::two_norm(u) ;
        for (typename vector_ortho_vec::iterator it = ortho_space_.begin(); it != ortho_space_.end(); it++)
            r -= ietl::dot((*it)[0], r) * (*it)[0] ;
        return r ;
    }
    // Routine doing deflation
    template <class Matrix, class VS>
    void davidson_modified<Matrix,VS>::update_orthospace(void)
    {
        for (size_t jcont = 0; jcont < n_root_found_; jcont++) {
            vector_type tmp   = vecspace_.return_orthovec(u_and_uA_[jcont][0], order_[n_root_found_], order_[jcont], site1_, site2_) ;
            vector_type tmpA  = vecspace_.return_orthovec(u_and_uA_[jcont][1], order_[n_root_found_], order_[jcont], site1_, site2_) ;
            vector_type tmpAA = vecspace_.return_orthovec(u_and_uA_[jcont][2], order_[n_root_found_], order_[jcont], site1_, site2_) ;
            for (size_t j = 0 ; j < ortho_space_.size() ; j++) {
                tmp   -= ietl::dot(ortho_space_[j][0], tmp) * ortho_space_[j][0] ;
                tmpA  -= ietl::dot(ortho_space_[j][0], tmp) * ortho_space_[j][1] ;
            }
            if (ietl::two_norm(tmp) > 1.0E-10) {
                tmp    /= ietl::two_norm(tmp);
                tmpA   /= ietl::two_norm(tmp);
                tmpAA  /= ietl::two_norm(tmp);
                std::vector< vector_type > junk ;
                junk.push_back(tmp)   ;
                junk.push_back(tmpA)  ;
                junk.push_back(tmpAA) ;
                ortho_space_.push_back(junk);
            }
        }
    }
    // Update the vector with the quantity to orthogonalize
    template <class Matrix, class VS>
    void davidson_modified<Matrix, VS>::update_u_and_uA(const vector_type &u, const vector_type &uA)
    {
        vector_type tmp1, tmp2, tmp3 ;
        tmp1 = u  / ietl::two_norm(u) ;
        tmp2 = uA / ietl::two_norm(u) ;
        tmp3 = apply_operator(tmp2) ;
        std::vector< vector_type > junk ;
        junk.push_back(tmp1) ;
        junk.push_back(tmp2) ;
        junk.push_back(tmp3) ;
        u_and_uA_.push_back(junk) ;
    }
    // Construction of the virtual function update_vspace
    template <class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::update_vspace(vector_set& V, vector_set& VA, vector_set& t)
    {
        size_t n_lin ;
        vector_set tA ;
        for (typename vector_ortho_vec::iterator it = ortho_space_.begin(); it != ortho_space_.end(); it++)
            t[0] -= ietl::dot((*it)[0], t[0]) * (*it)[0];
        tA[0] = apply_operator(t[0]) ;
        for (typename vector_ortho_vec::iterator it = ortho_space_.begin(); it != ortho_space_.end(); it++)
            tA[0] -= ietl::dot((*it)[0], tA[0]) * (*it)[0];
        n_lin = gram_schmidt_orthogonalizer_additional<vector_type,magnitude_type>(VA, V, tA, t);
    } ;
    // Definition of the virtual function precondition
    template<class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::precondition(vector_type &r, const vector_type &V, const vector_type& VA,
                                                     const magnitude_type &theta, const size_t& idx) {
        magnitude_type denom, x2, x1 ;
        x1 = ietl::dot(VA, r)/ietl::dot(VA,V) ;
        vector_type Vcpy = r - V * x1;
        bm_type &data = Vcpy.data();
        assert(shape_equal(data, Hdiag_[idx]));
        for (size_t b = 0; b < data.n_blocks(); ++b) {
            for (size_t i = 0; i < num_rows(data[b]); ++i) {
                for (size_t j = 0; j < num_cols(data[b]); ++j) {
                    denom = (omega_ - Hdiag_[idx][b](i, j)) - theta;
                    if (std::abs(denom))
                        data[b](i, j) /= denom;
                }
            }
        }
        x2 = ietl::dot(VA, Vcpy)/ietl::dot(VA,V);
        r = Vcpy - x2 * V ;
    } ;
    // Virtual function finalize_iteration
    template<class MATRIX, class VS>
    typename davidson_modified<MATRIX, VS>::vector_type davidson_modified<MATRIX, VS>::finalize_iteration
            (const vector_type &u, const vector_type &r, const size_t &n_restart,
             size_t &iter_dim, vector_set &V2, vector_set &VA)
    {
        vector_type result ;
        //if (iter_dim == n_restart){
        //    iter_dim = 0 ;
        //    V2.resize(0) ;
        //    VA.resize(0) ;
        //    result =  u ;
        //} else {
        //    result = r ;
        //}
        result = r ;
        return result  ;
    }
    // Routine to select the proper eigenpair
    template<class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::select_eigenpair(const vector_set& V2,
                                                         const vector_set& VA,
                                                         const matrix_numeric& Mevecs,
                                                         const vector_numeric& Mevals,
                                                         const size_t& dim,
                                                         vector_type& u,
                                                         vector_type& uA,
                                                         magnitude_type& theta)
    {
        // Initialization
        vector_type res ;
        // Main loop
        u = V2[0] * Mevecs(0,0) ;
        for (int i = 1; i < dim; ++i)
            u += V2[i] * Mevecs(i,0) ;
        uA = VA[0] * Mevecs(0,0) ;
        for (int i = 1; i < dim; ++i)
            uA += VA[i] * Mevecs(i,0) ;
        u  /= ietl::two_norm(uA) ;
        uA /= ietl::two_norm(uA) ;
        theta = Mevals[0] ;
        return ;
    }
    // Routine to print the header of the table
    template<class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::print_header_table(void) {
        printer_.print_header_table_simple() ;
    } ;
    //
    template<class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::print_endline(void) {
        printer_.print_endline_simple() ;
    } ;
    //
    template<class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::print_newline_table(const size_t& iter,
                                                            const size_t& size,
                                                            const magnitude_type& error,
                                                            const magnitude_type& energy)
    {
        printer_.print_newline_table_simple(iter, size, error, energy) ;
    } ;
}

#endif
