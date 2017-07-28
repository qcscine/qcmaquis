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
        typedef typename base::bm_type                bm_type;
        typedef typename base::magnitude_type         magnitude_type;
    	typedef typename base::matrix_numeric         matrix_numeric;
        typedef typename base::result_selection_type  result_selection_type ;
        typedef typename base::size_t                 size_t;
        typedef typename base::vector_set             vector_set;
        typedef typename base::vector_type            vector_type;
        using base::atol_ ;
        using base::Hdiag_ ;
        using base::matrix_ ;
        using base::n_sa_ ;
        using base::nsites_ ;
        using base::printer_ ;
        using base::site1_ ;
        using base::site2_ ;
        using base::vecspace_ ;
        // New constructors
        davidson_modified(const MATRIX &matrix, const VS &vec, const magnitude_type& omega, const int& nmin,
                          const int& nmax, const int& nsites, const int& site1, const int& site2)
                : base::davidson(matrix, vec, nmin, nmax, nsites, site1, site2) , omega_(omega) {};
        ~davidson_modified() {};
    private:
        // Private methods
        magnitude_type return_final(const magnitude_type &x) { return omega_-x; };
        result_selection_type select_eigenpair(const vector_set& V, const vector_set& VA, const matrix_numeric& eigvecs,
                                               const size_t& i, vector_set& u, vector_set& uA);
        vector_type apply_operator(const vector_type &x);
        vector_type finalize_iteration(const vector_type& u, const vector_type& r, const size_t& n_restart,
                                       size_t& iter_dim, vector_set& V2, vector_set& VA);
        void precondition(vector_type &r, const vector_type &V, const vector_type &VA, const magnitude_type &theta,
                          const size_t& idx);
        void update_vspace(vector_set &V, vector_set &VA, vector_set &t);
        // Printing-related methods
        void print_header_table(void) ;
        void print_endline(void) ;
        void print_newline_table(const size_t& iter, const size_t& size, const magnitude_type& error, const magnitude_type& energy,
                                 const magnitude_type& overlap=0.) ;
        void print_newline_table_energyonly(const magnitude_type& error, const magnitude_type& energy, const magnitude_type& overlap=0.) ;
        // Additional attributes
        magnitude_type omega_ ;
        vector_set V_additional_ ;
    };
    // Construction of the virtual function apply operator
    template <class MATRIX, class VS>
    typename davidson_modified<MATRIX, VS>::vector_type davidson_modified<MATRIX, VS>::apply_operator(const vector_type& x)
    {
        vector_type tmp ;
        ietl::mult(matrix_, x, tmp);
        tmp *= -1. ;
        tmp += omega_*x ;
        return tmp ;
    };
    // Construction of the virtual function update_vspace
    template <class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::update_vspace(vector_set& V, vector_set& VA, vector_set& t)
    {
        size_t n_lin_1 , n_lin_2 ;
        vector_set tA ;
        for (size_t i = 0 ; i < t.size() ; i++)
            tA.push_back(apply_operator(t[i])) ;
        //n_lin_1 = gram_schmidt_orthogonalizer_refinement<vector_type,magnitude_type>(V_additional_, t) ;
        n_lin_2 = gram_schmidt_orthogonalizer_additional<vector_type,magnitude_type>(VA, V, tA, t) ;
    } ;
    // Definition of the virtual function precondition
    template<class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::precondition(vector_type &r, const vector_type &V, const vector_type& VA, const magnitude_type &theta,
                                                     const size_t& idx) {
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
    typename davidson_modified<MATRIX, VS>::result_selection_type 
             davidson_modified<MATRIX, VS>::select_eigenpair(const vector_set& V2,
                                                             const vector_set& VA,
                                                             const matrix_numeric& Mevecs,
                                                             const size_t& dim,
                                                             vector_set& u,
                                                             vector_set& uA)
    {
        // Initialization
        assert(n_sa_ == u.size() && n_sa_ && uA.size()) ;
        result_selection_type res ;
        size_t n_eigen = std::min(dim, n_sa_) ;
        res.first = n_eigen ;
        // Main loop
        for (size_t k = 0 ; k < n_eigen ; k++) {
            u[k] = V2[0] * Mevecs(0,k);
            for (int i = 1; i < dim; ++i)
                u[k] += V2[i] * Mevecs(i,k);
            uA[k] = VA[0] * Mevecs(0,k);
            for (int i = 1; i < dim; ++i)
                uA[k] += VA[i] * Mevecs(i,k);
            uA[k] /= ietl::two_norm(u[k]);
            u[k] /= ietl::two_norm(u[k]);
            res.second.push_back(ietl::dot(uA[k],u[k])) ;
        }
        return res ;
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
                                                            const magnitude_type& energy,
                                                            const magnitude_type& overlap)
    {
        printer_.print_newline_table_simple(iter, size, error, energy) ;
    } ;
    //
    template<class MATRIX, class VS>
    void davidson_modified<MATRIX, VS>::print_newline_table_energyonly(const magnitude_type& error,
                                                                       const magnitude_type& energy,
                                                                       const magnitude_type& overlap )
    {
        printer_.print_newline_table_simple_onlyenergy(error, energy) ;
    } ;
}

#endif
