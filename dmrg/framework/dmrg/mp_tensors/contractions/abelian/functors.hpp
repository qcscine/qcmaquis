/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
        gemm_trim_left(A, B, C, ref_left_basis);
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
