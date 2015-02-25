/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef CONTRACTIONS_SU2_GEMM_HPP
#define CONTRACTIONS_SU2_GEMM_HPP

#include "dmrg/block_matrix/block_matrix.h"

namespace SU2 {

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm(block_matrix<Matrix1, SymmGroup> const & A,
              block_matrix<Matrix2, SymmGroup> const & B,
              block_matrix<Matrix3, SymmGroup> & C)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Index<SymmGroup>::const_iterator ci;

        C.clear();
        assert(B.basis().is_sorted());

        const_iterator BBbegin = B.basis().begin();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            std::pair<const_iterator, const_iterator>
              er = std::equal_range(BBbegin, B.basis().end(),
                dual_index_detail::QnBlock<SymmGroup>(A.basis().right_charge(k), SymmGroup::IdentityCharge, 0, 0), dual_index_detail::gt_row<SymmGroup>());

            for (const_iterator it = er.first; it != er.second; ++it)
            {
                std::size_t matched_block = std::distance(BBbegin, it);
                Matrix3 tmp(num_rows(A[k]), num_cols(B[matched_block]));
                gemm(A[k], B[matched_block], tmp);
                C.match_and_add_block(tmp, A.basis().left_charge(k), B.basis().right_charge(matched_block));
            }
        }
    }

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm_trim(block_matrix<Matrix1, SymmGroup> const & A,
                   block_matrix<Matrix2, SymmGroup> const & B,
                   block_matrix<Matrix3, SymmGroup> & C)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Index<SymmGroup>::const_iterator ci;

        C.clear();
        assert(B.basis().is_sorted());

        const_iterator BBbegin = B.basis().begin();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            std::pair<const_iterator, const_iterator>
              er = std::equal_range(BBbegin, B.basis().end(),
                dual_index_detail::QnBlock<SymmGroup>(A.basis().right_charge(k), SymmGroup::IdentityCharge, 0, 0), dual_index_detail::gt_row<SymmGroup>());

            for (const_iterator it = er.first; it != er.second; ++it)
            {
                std::size_t matched_block = std::distance(BBbegin, it);
                if (A.basis().left_charge(k) != B.basis().right_charge(matched_block)) continue;
                Matrix3 tmp(num_rows(A[k]), num_cols(B[matched_block]));
                gemm(A[k], B[matched_block], tmp);
                C.match_and_add_block(tmp, A.basis().left_charge(k), B.basis().right_charge(matched_block));
            }
        }
    }
}

#endif
