/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
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

    template <class T, class SymmGroup>
    T conjugate_correction(typename SymmGroup::charge lc, typename SymmGroup::charge rc, typename SymmGroup::subcharge tensor_spin)
    {
        assert( SymmGroup::spin(lc) >= 0);
        assert( SymmGroup::spin(rc) >= 0);

        typename SymmGroup::subcharge S = std::min(SymmGroup::spin(rc), SymmGroup::spin(lc));
        typename SymmGroup::subcharge spin_diff = SymmGroup::spin(rc) - SymmGroup::spin(lc);

        if (tensor_spin == 0)
        {
            return 1.;
        }
        else if (tensor_spin == 1)
        {
            if (spin_diff > 0)
                return -T( sqrt((S + 1.)/(S + 2.)) );

            else if (spin_diff < 0)
                return T( sqrt((S + 2.)/(S + 1.)) );
        }
        else if (tensor_spin == 2)
        {
            if (spin_diff > 0)
                return -T( sqrt( (S + 1.) / (S + 3.)) );

            else if (spin_diff < 0)
                return -T( sqrt((S + 3.) / (S + 1.)) );

            else
                return 1.;
        }
        else
            throw std::runtime_error("hermitian conjugate for reduced tensor operators only implemented up to rank 1");
    }

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm(block_matrix<Matrix1, SymmGroup> const & A,
              block_matrix<Matrix2, SymmGroup> const & B,
              block_matrix<Matrix3, SymmGroup> & C,
              int spin)
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
                std::size_t matched_block = std::distance(B_begin, it);
                if (!(spin == -1) && !::SU2::triangle(SymmGroup::spin(A.basis().left_charge(k)), spin, SymmGroup::spin(it->rc)))
                    continue;

                std::size_t c_block = C.find_block(A.basis().left_charge(k), it->rc);
                if (c_block == C.n_blocks())
                    c_block = C.insert_block(Matrix3(num_rows(A[k]), it->rs), A.basis().left_charge(k), it->rc);

                boost::numeric::bindings::blas::gemm(value_type(1), A[k], B[matched_block], value_type(1), C[c_block]);
            }
        }
    }

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm_trim_left(block_matrix<Matrix1, SymmGroup> const & A,
                        block_matrix<Matrix2, SymmGroup> const & B,
                        block_matrix<Matrix3, SymmGroup> & C,
                        std::vector<typename Matrix1::value_type> conj_scales = std::vector<typename Matrix1::value_type>())
    {
        typedef typename SymmGroup::charge charge;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Matrix3::value_type value_type;

        if (conj_scales.size() != A.n_blocks())
            conj_scales = std::vector<value_type>(A.n_blocks(), 1.);

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

                boost::numeric::bindings::blas::gemm(conj_scales[k], A[k], B[matched_block], value_type(1), C[c_block]);
            }
        }
    }

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm_trim_right(block_matrix<Matrix1, SymmGroup> const & A,
                         block_matrix<Matrix2, SymmGroup> const & B,
                         block_matrix<Matrix3, SymmGroup> & C,
                         std::vector<typename Matrix1::value_type> conj_scales = std::vector<typename Matrix1::value_type>())
    {
        typedef typename SymmGroup::charge charge;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Matrix3::value_type value_type;

        if (conj_scales.size() != B.n_blocks())
            conj_scales = std::vector<value_type>(B.n_blocks(), 1.);

        C.clear();
        assert(B.basis().is_sorted());

        const_iterator B_begin = B.basis().begin();
        const_iterator B_end = B.basis().end();
        Index<SymmGroup> A_right_basis = A.right_basis();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            charge ar = A.basis().right_charge(k);

            const_iterator it = B.basis().left_lower_bound(ar);

            for ( ; it != B_end && it->lc == ar; ++it)
            {
                std::size_t matched_block = std::distance(B_begin, it);
                if (!A_right_basis.has(it->rc)) continue;

                std::size_t c_block = C.find_block(A.basis().left_charge(k), it->rc);
                if (c_block == C.n_blocks())
                    c_block = C.insert_block(Matrix3(num_rows(A[k]), it->rs), A.basis().left_charge(k), it->rc);

                boost::numeric::bindings::blas::gemm(conj_scales[matched_block], A[k], B[matched_block], value_type(1), C[c_block]);
            }
        }
    }

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm_trim(block_matrix<Matrix1, SymmGroup> const & A,
                   block_matrix<Matrix2, SymmGroup> const & B,
                   block_matrix<Matrix3, SymmGroup> & C,
                   std::vector<typename Matrix1::value_type> conj_scales,
                   bool conjugate_a)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Matrix3::value_type value_type;

        assert(B.basis().is_sorted());
        assert( (conjugate_a && A.n_blocks() == conj_scales.size()) || (!conjugate_a && B.n_blocks() == conj_scales.size()));
        C.clear();

        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            charge al = A.basis().left_charge(k);
            charge ar = A.basis().right_charge(k);

            std::size_t matched_block = B.basis().position(ar, al);
            if (matched_block == B.n_blocks()) continue;

            std::size_t c_block = C.find_block(al, al);
            if (c_block == C.n_blocks())
                c_block = C.insert_block(Matrix3(num_rows(A[k]), num_cols(B[matched_block])), al, al);

            //boost::numeric::bindings::blas::gemm(value_type(1), A[k], B[matched_block], value_type(1), C[c_block]);
            boost::numeric::bindings::blas::gemm(conj_scales[ (conjugate_a) ? k : matched_block], A[k], B[matched_block],
                                                 value_type(1), C[c_block]);
        }
    }
}

#endif
