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
    void gemm_trim_left(block_matrix<Matrix1, SymmGroup> const & A,
                        block_matrix<Matrix2, SymmGroup> const & B,
                        block_matrix<Matrix3, SymmGroup> & C)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Matrix3::value_type value_type;

        C.clear();
        assert(B.basis().is_sorted());

        const_iterator B_begin = B.basis().begin();
        const_iterator B_end = B.basis().end();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            if (!B.basis().left_has(A.basis().left_charge(k))) continue;

            charge ar = A.basis().right_charge(k);
            const_iterator it = B.basis().left_lower_bound(ar);

            for ( ; it != B_end && it->lc == ar; ++it)
            {
                std::size_t matched_block = std::distance(B_begin, it);

                std::size_t c_block = C.find_block(A.basis().left_charge(k), it->rc);
                if (c_block == C.n_blocks())
                    c_block = C.insert_block(Matrix3(num_rows(A[k]), it->rs), A.basis().left_charge(k), it->rc);

                boost::numeric::bindings::blas::gemm(value_type(1), A[k], B[matched_block], value_type(1), C[c_block]);
            }

            //for ( ; it != B_end && it->lc == ar; ++it)
            //{
            //    std::size_t matched_block = std::distance(B_begin, it);
            //    Matrix3 tmp(num_rows(A[k]), it->rs);
            //    gemm(A[k], B[matched_block], tmp);
            //    C.match_and_add_block(tmp, A.basis().left_charge(k), it->rc);
            //}
        }
    }

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm_trim(block_matrix<Matrix1, SymmGroup> const & A,
                   block_matrix<Matrix2, SymmGroup> const & B,
                   block_matrix<Matrix3, SymmGroup> & C)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Matrix3::value_type value_type;

        C.clear();
        assert(B.basis().is_sorted());

        const_iterator B_begin = B.basis().begin();
        const_iterator B_end = B.basis().end();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            charge ar = A.basis().right_charge(k);
            const_iterator it = B.basis().left_lower_bound(ar);

            for ( ; it != B_end && it->lc == ar; ++it)
            {
                if (A.basis().left_charge(k) != it->rc) continue;

                std::size_t matched_block = std::distance(B_begin, it);

                std::size_t c_block = C.find_block(A.basis().left_charge(k), it->rc);
                if (c_block == C.n_blocks())
                    c_block = C.insert_block(Matrix3(num_rows(A[k]),it->rs), A.basis().left_charge(k), it->rc);

                boost::numeric::bindings::blas::gemm(value_type(1), A[k], B[matched_block], value_type(1), C[c_block]);
            }

            //for ( ; it != B_end && it->lc == ar; ++it)
            //{
            //    if (A.basis().left_charge(k) != it->rc) continue;
            //    std::size_t matched_block = std::distance(B_begin, it);
            //    Matrix3 tmp(num_rows(A[k]), it->rs);
            //    gemm(A[k], B[matched_block], tmp);
            //    C.match_and_add_block(tmp, A.basis().left_charge(k), it->rc);
            //}
        }
    }
}

#endif
