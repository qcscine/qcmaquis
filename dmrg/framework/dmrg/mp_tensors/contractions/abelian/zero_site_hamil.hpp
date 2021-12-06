/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Laboratory for Physical Chemistry, ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef ABELIAN_ZERO_SITE_HAMIL
#define ABELIAN_ZERO_SITE_HAMIL

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"

namespace contraction {

template<class Matrix, class OtherMatrix, class SymmGroup, class SymmType>
block_matrix<Matrix, SymmGroup>
Engine<Matrix, OtherMatrix, SymmGroup, SymmType>::zerosite_hamil2(block_matrix<Matrix, SymmGroup> ket_tensor,
                                                                  Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                                                                  MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right)
{
    return zerosite_hamil2_kernel(ket_tensor, left, right, mpo_left, mpo_right);
}

template<class Matrix, class OtherMatrix, class SymmGroup, class SymmType>
block_matrix<Matrix, SymmGroup>
Engine<Matrix, OtherMatrix, SymmGroup, SymmType>::zerosite_hamil2_kernel(block_matrix<Matrix, SymmGroup> ket_bm,
                                                                         Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                                                                         MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right)
{
    using charge = typename SymmGroup::charge;
    using index_type = typename MPOTensor<Matrix, SymmGroup>::index_type;
    BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, abelian::Gemms> t(ket_bm, left, mpo_right, false);
    block_matrix<Matrix, SymmGroup> ret;
    index_type loop_max = mpo_right.row_dim();
    omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
        block_matrix<Matrix, SymmGroup> tmp, local;
        if (mpo_right.herm_info.left_skip(b2)) {
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
