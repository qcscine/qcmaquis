/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
 *               2016-2016 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2018 by Alberto Baiardi <abaiardi@ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

#ifndef CONTRACTIONS_ABELIAN_DETAIL_HPP
#define CONTRACTIONS_ABELIAN_DETAIL_HPP

#include "dmrg/block_matrix/block_matrix.h"

namespace contraction {
namespace abelian {
namespace detail {

    template<class Matrix1, class SymmGroup>
    DualIndex<SymmGroup> gemm_trim_left_basis(block_matrix<Matrix1, SymmGroup> const & A,
                                              DualIndex<SymmGroup> const & B_basis,
                                              Index<SymmGroup> const & ref_left_basis)
    {
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;

        DualIndex<SymmGroup> ret;

        const_iterator B_begin = B_basis.begin();
        const_iterator B_end = B_basis.end();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            // if (!B_basis.left_has(A.basis().left_charge(k))) continue;
            if (!ref_left_basis.has(A.basis().left_charge(k))) continue;

            typename SymmGroup::charge ar = A.basis().right_charge(k);
            const_iterator it = B_basis.left_lower_bound(ar);

            for ( ; it != B_end && it->lc == ar; ++it)
                if (!ret.has(A.basis().left_charge(k), it->rc))
                    ret.insert(typename DualIndex<SymmGroup>::value_type(A.basis().left_charge(k), it->rc, A.basis().left_size(k), it->rs));
        }
        return ret;
    }
    template<class Matrix2, class SymmGroup>
    DualIndex<SymmGroup> gemm_trim_right_basis(DualIndex<SymmGroup> const & A_basis,
                                               block_matrix<Matrix2, SymmGroup> const & B,
                                               Index<SymmGroup> const & ref_right_basis)
    {
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;

        DualIndex<SymmGroup> ret;

        const_iterator B_begin = B.basis().begin();
        const_iterator B_end = B.basis().end();
        Index<SymmGroup> A_right_basis(A_basis.size());
        for (size_t k = 0; k < A_basis.size(); ++k)
            A_right_basis[k] = std::make_pair(A_basis.right_charge(k), A_basis.right_size(k));

        for (std::size_t k = 0; k < A_basis.size(); ++k) {

            typename SymmGroup::charge ar = A_basis.right_charge(k);
            const_iterator it = B.basis().left_lower_bound(ar);

            for ( ; it != B_end && it->lc == ar; ++it)
            {
                // if (!A_right_basis.has(it->rc)) continue; // trim
                if (!ref_right_basis.has(it->rc)) continue; // trim
                if (!ret.has(A_basis.left_charge(k), it->rc))
                    ret.insert(typename DualIndex<SymmGroup>::value_type(A_basis.left_charge(k), it->rc, A_basis.left_size(k), it->rs));
            }
        }
        return ret;
    }



    template<class Matrix, class OtherMatrix, class SymmGroup>
    DualIndex<SymmGroup> T_basis_left(Boundary<OtherMatrix, SymmGroup> const & boundary,
                                      common::BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemms> const & mult_mps,
                                      MPOTensor<Matrix, SymmGroup> const & mpo,
                                      DualIndex<SymmGroup> const & mps_basis,
                                      Index<SymmGroup> const & ref_left_basis,
                                      typename MPOTensor<Matrix, SymmGroup>::index_type b)
    {
        if (mpo.num_row_non_zeros(b) == 1)
            if (mpo.herm_info.left_skip(b))
                return gemm_trim_left_basis(conjugate(boundary[mpo.herm_info.left_conj(b)]), mps_basis, ref_left_basis);
            else
                return gemm_trim_left_basis(transpose(boundary[b]), mps_basis, ref_left_basis);
        else
            return mult_mps[b].basis();
    }
    template<class Matrix, class OtherMatrix, class SymmGroup>
    DualIndex<SymmGroup> T_basis_right(Boundary<OtherMatrix, SymmGroup> const & boundary,
                                       common::MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, Gemms> const & mult_mps,
                                       MPOTensor<Matrix, SymmGroup> const & mpo,
                                       DualIndex<SymmGroup> const & mps_basis,
                                       Index<SymmGroup> const & ref_right_basis,
                                       typename MPOTensor<Matrix, SymmGroup>::index_type b)
    {
        if (mpo.num_col_non_zeros(b) == 1)
            if (mpo.herm_info.right_skip(b))
                return gemm_trim_right_basis(mps_basis, transpose(conjugate(boundary[mpo.herm_info.right_conj(b)])), ref_right_basis);
            else
                return gemm_trim_right_basis(mps_basis, boundary[b], ref_right_basis);
        else
            return mult_mps[b].basis();
    }

}
} // namespace abelian
} // namespace contraction

#endif
