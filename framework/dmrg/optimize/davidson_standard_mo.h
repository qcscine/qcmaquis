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
        typedef typename base::result_selection_type result_selection_type ;
        typedef typename base::size_t                size_t;
        typedef typename base::vector_set            vector_set;
        typedef typename base::vector_type           vector_type;
        using base::Hdiag_;
        using base::matrix_;
        using base::nsites_ ;
        using base::n_sa_ ;
        using base::site1_;
        using base::site2_;
        using base::printer_ ;
        using base::v_guess_ ;
        // New constructors
        davidson_standard_mo(const MATRIX &matrix, const VS &vec, const pov_vec poverlap,
                             const int& nmin, const int& nmax, const int& nsites, const int& site1,
                             const int& site2, const int& root_homing_type)
                : base::davidson(matrix, vec, nmin, nmax, nsites, site1, site2),
                  pov_(poverlap), root_homing_type_(root_homing_type) {};
        ~davidson_standard_mo() {};
    private:
        // Private methods
        void precondition(vector_type &r, const vector_type &V, const vector_type& VA, const magnitude_type &theta);
        result_selection_type select_eigenpair(const vector_set& V, const vector_set& VA, const matrix_numeric& eigvecs,
                                               const size_t& i, vector_set& u, vector_set& uA);
        void update_vspace(vector_set &V, vector_set &VA, vector_set &t, const size_t& dim);
        vector_type apply_operator(const vector_type &x);
        magnitude_type return_final(const magnitude_type &x) { return x; };
        vector_type finalize_iteration(const vector_type& u, const vector_type& r, const size_t& n_restart,
                                       size_t& iter_dim, vector_set& V2, vector_set& VA);
        double compute_overlap(const vector_type& vec_test, const size_t& idx) ;
        // Private attribute
        int root_homing_type_ ;
        pov_vec  pov_ ;
    };
    // Definition of the virtual function update_vspace
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::update_vspace(vector_set &V, vector_set &VA,
                                                                                 vector_set &t, const size_t& dim)
    {
        for (size_t j = 0 ; j < n_sa_ ; j++) {
            for (size_t i = 0; i < dim; i++)
                t[j] -= ietl::dot(V[i], t[j]) * V[i];
            V.push_back(t / ietl::two_norm(t[j]));
            VA.resize(V.size());
            VA[V.size()-1] = apply_operator(V[V.size()-1]) ;
        }
    };
    // Definition of the virtual function apply_operator
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    typename davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::vector_type
             davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::apply_operator(const vector_type &x) {
        vector_type tmp;
        ietl::mult(matrix_, x, tmp);
        return tmp;
    };
    // Definition of the virtual function precondition
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::precondition(vector_type &r, const vector_type &V, const vector_type& VA, const magnitude_type &theta) {
        magnitude_type denom, x2, x1 = ietl::dot(V, r);
        vector_type Vcpy = r - V * x1;
        bm_type &data = Vcpy.data();
        assert(shape_equal(data, Hdiag_));
        for (size_t b = 0; b < data.n_blocks(); ++b) {
            for (size_t i = 0; i < num_rows(data[b]); ++i) {
                for (size_t j = 0; j < num_cols(data[b]); ++j) {
                    denom = Hdiag_[b](i, j) - theta;
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
        if (iter_dim == n_restart){
            iter_dim = 0 ;
            V2.resize(0) ;
            VA.resize(0) ;
            result =  u ;
        } else {
            result = r ;
        }
        return result  ;
    }
    // Routine to select the proper eigenpair
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    typename davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::result_selection_type 
             davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::select_eigenpair(const vector_set& V, 
                                                                                        const vector_set& VA, 
                                                                                        const matrix_numeric& Mevecs,
                                                                                        const size_t& dim, 
                                                                                        vector_set& u, 
                                                                                        vector_set& uA)
    {
        int idx = 0 ;
        double scr ;
        size_t n_eigen = std::min(n_sa_,dim) ;
        std::vector<double> overlaps(dim);
        result_selection_type res ;
        vector_type u_local ;
        res.first = n_eigen ;
        //
        for (size_t k = 0 ; k < n_eigen ; k++) {
            for (size_t i = 0; i < dim; ++i) {
                // Conversion to the original basis
                u_local = Mevecs(0,i) * V[0];
                for (size_t j = 1; j < dim; ++j)
                    u_local += Mevecs(j,i) * V[j];
                overlaps[i] = compute_overlap(u_local,k) ;
            }
            for (size_t i = 1; i < dim; ++i) {
                if (overlaps[i] > overlaps[idx])
                    idx = i;
            }
            //
            std::cout << "Overlap - " << overlaps[idx] << std::endl ;
            res.second.push_back(overlaps[idx]) ;
            //
            u[k]  = V[0]*Mevecs(0,idx) ;
            uA[k] = VA[0]*Mevecs(0,idx);
            for (int i = 1; i < dim; ++i) {
                u[k]  += Mevecs(i, idx) * V[i];
                uA[k] += Mevecs(i, idx) * VA[i];
            }
            uA[k] /= ietl::two_norm(u[k]) ;
            u[k]  /= ietl::two_norm(u[k]) ;
        }
        return res ;
    }
    // Routine to compute the overlaps
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    double davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::compute_overlap(const vector_type &vec_test, const size_t& idx)
    {
        double ret, scr ;
        if (root_homing_type_ == 1)
            if (nsites_ == 1)
                ret = pov_[idx].overlap(vec_test / ietl::two_norm(vec_test), site1_);
            else
                ret = pov_[idx].overlap(vec_test / ietl::two_norm(vec_test), site1_, site2_);
        else
            ret = ietl::dot(vec_test, v_guess_[0]) / ietl::two_norm(vec_test);
        return fabs(ret) ;
    }
}

#endif
