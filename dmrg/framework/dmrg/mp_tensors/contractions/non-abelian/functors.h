/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_SU2_FUNCTORS_HPP
#define CONTRACTIONS_SU2_FUNCTORS_HPP


#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/contractions/non-abelian/gemm.hpp"

namespace SU2 {

    struct su2gemm
    {
        template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
        void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                        block_matrix<Matrix2, SymmGroup> const & B,
                        block_matrix<Matrix3, SymmGroup> & C,
                        int spin = -1)
        {
            SU2::gemm(A,B,C, spin);
        }
    };

    struct su2gemm_trim_left
    {
        template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
        void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                        block_matrix<Matrix2, SymmGroup> const & B,
                        block_matrix<Matrix3, SymmGroup> & C,
                        std::vector<typename Matrix1::value_type> scales = std::vector<typename Matrix1::value_type>())
        {
            SU2::gemm_trim_left(A, B, C, scales);
        }

        template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
        void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                        block_matrix<Matrix2, SymmGroup> const & B,
                        block_matrix<Matrix3, SymmGroup> & C,
                        Index<SymmGroup> const & ref_left_basis,
                        std::vector<typename Matrix1::value_type> scales = std::vector<typename Matrix1::value_type>())
        {
            SU2::gemm_trim_left(A, B, C, ref_left_basis, scales);
        }
    };

    struct su2gemm_trim_right
    {
        template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
        void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                        block_matrix<Matrix2, SymmGroup> const & B,
                        block_matrix<Matrix3, SymmGroup> & C,
                        std::vector<typename Matrix1::value_type> scales = std::vector<typename Matrix1::value_type>())
        {
            SU2::gemm_trim_right(A,B,C, scales);
        }

        template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
        void operator()(block_matrix<Matrix1, SymmGroup> const & A,
                        block_matrix<Matrix2, SymmGroup> const & B,
                        block_matrix<Matrix3, SymmGroup> & C,
                        Index<SymmGroup> const & ref_right_basis,
                        std::vector<typename Matrix1::value_type> scales = std::vector<typename Matrix1::value_type>())
        {
            SU2::gemm_trim_right(A, B, C, ref_right_basis, scales);
        }
    };

    struct SU2Gemms
    {
        typedef su2gemm gemm;
        typedef su2gemm_trim_left  gemm_trim_left;
        typedef su2gemm_trim_right gemm_trim_right;
    };

} // namespace SU2

#endif
