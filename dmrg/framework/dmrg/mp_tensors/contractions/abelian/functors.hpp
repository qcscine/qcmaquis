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

#ifndef CONTRACTIONS_ABELIAN_FUNCTORS_HPP
#define CONTRACTIONS_ABELIAN_FUNCTORS_HPP


#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

namespace contraction {
namespace abelian {

        struct gemm_functor
        {
            template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
            void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                            block_matrix<Matrix2, SymmGroup> const & B,
                            block_matrix<Matrix3, SymmGroup> & C,
                            int spin = -1)
            {
                gemm(A,B,C);
            }
        };

        struct gemm_trim_left_functor
        {
            template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
            void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                            block_matrix<Matrix2, SymmGroup> const & B,
                            block_matrix<Matrix3, SymmGroup> & C,
                            std::vector<typename Matrix1::value_type> scales = std::vector<typename Matrix1::value_type>())
            {
                gemm_trim_left(A,B,C, B.left_basis());
            }

            template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
            void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                            block_matrix<Matrix2, SymmGroup> const & B,
                            block_matrix<Matrix3, SymmGroup> & C,
                            Index<SymmGroup> const & ref_left_basis,
                            std::vector<typename Matrix1::value_type> scales = std::vector<typename Matrix1::value_type>())
            {
                gemm_trim_left(A,B,C, ref_left_basis);
            }

        };

        struct gemm_trim_right_functor
        {
            template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
            void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                            block_matrix<Matrix2, SymmGroup> const & B,
                            block_matrix<Matrix3, SymmGroup> & C,
                            std::vector<typename Matrix1::value_type> scales = std::vector<typename Matrix1::value_type>())
            {
                gemm_trim_right(A,B,C, A.right_basis());
            }

            template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
            void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                            block_matrix<Matrix2, SymmGroup> const & B,
                            block_matrix<Matrix3, SymmGroup> & C,
                            Index<SymmGroup> const & ref_right_basis,
                            std::vector<typename Matrix1::value_type> scales = std::vector<typename Matrix1::value_type>())
            {
                gemm_trim_right(A,B,C, ref_right_basis);
            }
        };

        struct Gemms
        {
            typedef gemm_functor gemm;
            typedef gemm_trim_left_functor gemm_trim_left;
            typedef gemm_trim_right_functor gemm_trim_right;
        };

} // namespace abelian
} // namespace contraction

#endif
