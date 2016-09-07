/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
 *               2016-2016 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

    template<class Matrix1, class Matrix2, class SymmGroup>
    DualIndex<SymmGroup> gemm_trim_left_basis(block_matrix<Matrix1, SymmGroup> const & A,
                                              block_matrix<Matrix2, SymmGroup> const & B)
    {
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;

        DualIndex<SymmGroup> ret;

        const_iterator B_begin = B.basis().begin();
        const_iterator B_end = B.basis().end();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            if (!B.basis().left_has(A.basis().left_charge(k))) continue;

            typename SymmGroup::charge ar = A.basis().right_charge(k);
            const_iterator it = B.basis().left_lower_bound(ar);

            for ( ; it != B_end && it->lc == ar; ++it)
                if (!ret.has(A.basis().left_charge(k), it->rc))
                    ret.insert(typename DualIndex<SymmGroup>::value_type(A.basis().left_charge(k), it->rc, A.basis().left_size(k), it->rs));
        }
        return ret;
    }
    template<class Matrix1, class Matrix2, class SymmGroup>
    DualIndex<SymmGroup> gemm_trim_right_basis(block_matrix<Matrix1, SymmGroup> const & A,
                                               block_matrix<Matrix2, SymmGroup> const & B)
    {
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;

        DualIndex<SymmGroup> ret;

        const_iterator B_begin = B.basis().begin();
        const_iterator B_end = B.basis().end();
        Index<SymmGroup> A_right_basis = A.right_basis();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            typename SymmGroup::charge ar = A.basis().right_charge(k);
            const_iterator it = B.basis().left_lower_bound(ar);

            for ( ; it != B_end && it->lc == ar; ++it)
            {
                if (!A_right_basis.has(it->rc)) continue; // trim
                if (!ret.has(A.basis().left_charge(k), it->rc))
                    ret.insert(typename DualIndex<SymmGroup>::value_type(A.basis().left_charge(k), it->rc, A.basis().left_size(k), it->rs));
            }
        }
        return ret;
    }



    template<class Matrix, class OtherMatrix, class SymmGroup>
    DualIndex<SymmGroup> T_basis_left(Boundary<OtherMatrix, SymmGroup> const & boundary,
                                      common::BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemms> const & mult_mps,
                                      MPOTensor<Matrix, SymmGroup> const & mpo,
                                      MPSTensor<Matrix, SymmGroup> const & mps,
                                      typename MPOTensor<Matrix, SymmGroup>::index_type b)
    {
        if (mpo.num_row_non_zeros(b) == 1)
            return gemm_trim_left_basis(transpose(boundary[b]), mps.data());
        else
            return mult_mps[b].basis();
    }
    template<class Matrix, class OtherMatrix, class SymmGroup>
    DualIndex<SymmGroup> T_basis_right(Boundary<OtherMatrix, SymmGroup> const & boundary,
                                       common::MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, Gemms> const & mult_mps,
                                       MPOTensor<Matrix, SymmGroup> const & mpo,
                                       MPSTensor<Matrix, SymmGroup> const & mps,
                                       typename MPOTensor<Matrix, SymmGroup>::index_type b)
    {
        if (mpo.num_col_non_zeros(b) == 1)
            return gemm_trim_right_basis(mps.data(), boundary[b]);
        else
            return mult_mps[b].basis();
    }



    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup> const * T_left(Boundary<OtherMatrix, SymmGroup> const & boundary,
                                                   common::BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemms> const & mult_mps,
                                                   block_matrix<Matrix, SymmGroup> & local,
                                                   MPOTensor<Matrix, SymmGroup> const & mpo,
                                                   MPSTensor<Matrix, SymmGroup> const & mps,
                                                   typename MPOTensor<Matrix, SymmGroup>::index_type b)
    {
        block_matrix<Matrix, SymmGroup> const * Tp = &local;
        if (mpo.num_row_non_zeros(b) == 1)
            gemm_trim_left(transpose(boundary[b]), mps.data(), local);
        else
            Tp = &mult_mps[b];

        return Tp;
    }
    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup> const * T_right(Boundary<OtherMatrix, SymmGroup> const & boundary,
                                                    common::MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, Gemms> const & mult_mps,
                                                    block_matrix<Matrix, SymmGroup> & local,
                                                    MPOTensor<Matrix, SymmGroup> const & mpo,
                                                    MPSTensor<Matrix, SymmGroup> const & mps,
                                                    typename MPOTensor<Matrix, SymmGroup>::index_type b)
    {
        block_matrix<Matrix, SymmGroup> const * Tp = &local;
        if (mpo.num_col_non_zeros(b) == 1)
            gemm_trim_right(mps.data(), boundary[b], local);
        else
            Tp = &mult_mps[b];

        return Tp;
    }
}
} // namespace abelian
} // namespace contraction

#endif
