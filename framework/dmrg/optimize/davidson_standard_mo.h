/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2015 by Alberto Baiardi <alberto.baiardi@sns.it>
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
        typedef davidson<MATRIX, VS> base;
        typedef typename base::bm_type        bm_type;
        typedef typename base::magnitude_type magnitude_type;
        typedef typename base::matrix_numeric matrix_numeric;
        typedef typename base::vector_set     vector_set;
        typedef typename base::vector_type    vector_type;
        typedef typename base::size_t         size_t;
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap ;
        using base::matrix_;
        using base::vecspace_;
        using base::atol_;
        using base::Hdiag_;
        using base::site_;
        // New constructors
        davidson_standard_mo(const MATRIX &matrix, const VS &vec, const int& site, const partial_overlap poverlap,
                             const int& nmin, const int& nmax)
                : base::davidson(matrix, vec, site, nmin, nmax), pov_(poverlap) {};
        ~davidson_standard_mo() {};
    private:
        // Private methods
        void precondition(vector_type &r, const vector_type &V, const magnitude_type &theta);
        void select_eigenpair(const vector_set& V, const vector_set& VA, const matrix_numeric& eigvecs,
                              const size_t& i, vector_type& u, vector_type& uA);
        void update_vspace(vector_set &V, vector_set &VA, vector_type &t, std::size_t dim);
        vector_type apply_operator(const vector_type &x);
        magnitude_type return_final(const magnitude_type &x) { return x; };
        vector_type finalize_iteration(const vector_type& u, const vector_type& r, const size_t& n_restart,
                                       size_t& iter_dim, vector_set& V2, vector_set& VA);
        // Private attribute
        partial_overlap  pov_ ;
    };
    // Definition of the virtual function update_vspace
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::update_vspace(vector_set &V, vector_set &VA, vector_type &t, std::size_t dim) {
        magnitude_type tau = ietl::two_norm(t);
        for (int i = 0; i < dim; i++)
            t -= ietl::dot(V[i], t) * V[i];
        V.push_back(t / ietl::two_norm(t));
        VA.resize(V.size());
        ietl::mult(matrix_, V[V.size() - 1], VA[V.size() - 1]);
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
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::precondition(vector_type &r, const vector_type &V, const magnitude_type &theta) {
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
            (const vector_type &u, const vector_type &r, const size_t &n_restart, size_t &iter_dim, vector_set &V2, vector_set &VA)
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
    void davidson_standard_mo<MATRIX, VS, OtherMatrix, SymmGroup>::select_eigenpair(const vector_set& V, const vector_set& VA, const matrix_numeric& Mevecs,
                                                                                    const size_t& dim, vector_type& u, vector_type& uA)
    {
        int idx = 0 ;
        vector_type u_local ;
        std::vector<double> overlaps ;
        overlaps.resize(dim) ;
        //
        for (int i = 0; i < dim; ++i) {
            // Conversion to the original basis
            u_local = Mevecs(0,i) * V[0];
            for (int j = 1; j < dim; ++j)
                u_local += Mevecs(j,i) * V[j];
            double scr = pov_.overlap(u_local, site_);
            overlaps[i] = fabs(scr);
        }
        for (int i = 1; i < dim; ++i) {
            if (overlaps[i] > overlaps[idx])
                idx = i;
        }
        //
        std::cout << "Overlap - " << overlaps[idx] << std::endl ;
        //
        u  = V[0]*Mevecs(0,idx) ;
        uA = VA[0]*Mevecs(0,idx);
        for (int i = 1; i < dim; ++i) {
            u  += Mevecs(i, idx) * V[i];
            uA += Mevecs(i, idx) * VA[i];
        }
        uA /= ietl::two_norm(u) ;
        u  /= ietl::two_norm(u) ;
    }
}

#endif
