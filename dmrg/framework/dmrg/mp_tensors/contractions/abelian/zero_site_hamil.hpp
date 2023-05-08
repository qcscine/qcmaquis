/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef ABELIAN_ZERO_SITE_HAMIL
#define ABELIAN_ZERO_SITE_HAMIL

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"

namespace contraction {

template<class Matrix, class OtherMatrix, class SymmGroup, class SymmType>
block_matrix<Matrix, SymmGroup>
Engine<Matrix, OtherMatrix, SymmGroup, SymmType>::zerosite_hamil2(block_matrix<Matrix, SymmGroup> ket_tensor,
                                                                  Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                                                                  MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right,
                                                                  bool isHermitian)
{
    return zerosite_hamil2_kernel(ket_tensor, left, right, mpo_left, mpo_right, isHermitian);
}

template<class Matrix, class OtherMatrix, class SymmGroup, class SymmType>
block_matrix<Matrix, SymmGroup>
Engine<Matrix, OtherMatrix, SymmGroup, SymmType>::zerosite_hamil2_kernel(block_matrix<Matrix, SymmGroup> ket_bm,
                                                                         Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                                                                         MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right,
                                                                         bool isHermitian)
{
    using charge = typename SymmGroup::charge;
    using index_type = typename MPOTensor<Matrix, SymmGroup>::index_type;
    BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, abelian::Gemms> t(ket_bm, left, mpo_right, isHermitian, false);
    block_matrix<Matrix, SymmGroup> ret;
    index_type loop_max = mpo_right.row_dim();
    omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
        block_matrix<Matrix, SymmGroup> tmp, local;
        if (mpo_right.herm_info.left_skip(b2) && isHermitian) {
            gemm(t.at(b2, local), adjoint(right[mpo_left.herm_info.right_conj(b2)]), tmp);
        }
        else {
            gemm(t.at(b2, local), right[b2], tmp);
        }
    parallel_critical
        for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
            ret.match_and_add_block(tmp[k], tmp.basis().left_charge(k), tmp.basis().right_charge(k));
    });
    return ret;
}

} // namespace contraction

#endif
