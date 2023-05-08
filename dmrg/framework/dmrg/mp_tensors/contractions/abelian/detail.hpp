/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_ABELIAN_DETAIL_HPP
#define CONTRACTIONS_ABELIAN_DETAIL_HPP

#include "dmrg/block_matrix/block_matrix.h"

namespace contraction {
namespace abelian {
namespace detail {

/**
 * @brief Returns the indices associated with the product of two [block_matrix] objects.
 * 
 * Same as [gemm_trim_left], but does not do any "real" calculation, just returns 
 * the DualIndex object associated with the product of two block_matrix objects.
 * Note that the second matrix is not passed as input, but is passed only as a
 * DualIndex object.
 * 
 * @param A First Matrix
 * @param B_basis DualIndex associated with the second matrix.
 * @param refBasis reference DualIndex used for trimming (as in [gemm_trim_left])
 * @return DualIndex<SymmGroup> Output DualIndex associated with A*B
 */
template<class Matrix1, class SymmGroup>
DualIndex<SymmGroup> gemm_trim_left_basis(block_matrix<Matrix1, SymmGroup> const & A,
                                          DualIndex<SymmGroup> const & B_basis,
                                          DualIndex<SymmGroup> const & refBasis)
{
    typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
    DualIndex<SymmGroup> ret;
    const_iterator B_begin = B_basis.begin();
    const_iterator B_end = B_basis.end();
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        if (!refBasis.left_has(A.basis().left_charge(k))) 
            continue;
        typename SymmGroup::charge ar = A.basis().right_charge(k);
        const_iterator it = B_basis.left_lower_bound(ar);
        for ( ; it != B_end && it->lc == ar; ++it)
            if (!ret.has(A.basis().left_charge(k), it->rc))
                ret.insert(typename DualIndex<SymmGroup>::value_type(A.basis().left_charge(k), it->rc, A.basis().left_size(k), it->rs));
    }
    return ret;
}


/**
 * @brief Returns the indices associated with the product of two [block_matrix] objects.
 * 
 * This is the right counterpart of [gemm_trim_left_basis]
 * 
 * @param A First Matrix
 * @param B_basis DualIndex associated with the second matrix.
 * @param refBasis reference DualIndex used for right trimming (as in [gemm_trim_right])
 * @return DualIndex<SymmGroup> Output DualIndex associated with A*B
 */
template<class Matrix2, class SymmGroup>
DualIndex<SymmGroup> gemm_trim_right_basis(DualIndex<SymmGroup> const & A_basis,
                                           block_matrix<Matrix2, SymmGroup> const & B,
                                           DualIndex<SymmGroup> const & refBasis)
{
    typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
    DualIndex<SymmGroup> ret;
    const_iterator B_begin = B.basis().begin();
    const_iterator B_end = B.basis().end();
    Index<SymmGroup> A_right_basis(refBasis.size());
    for (size_t k = 0; k < refBasis.size(); ++k) 
        A_right_basis[k] = std::make_pair(refBasis.right_charge(k), refBasis.right_size(k));
    for (std::size_t k = 0; k < A_basis.size(); ++k) {
        typename SymmGroup::charge ar = A_basis.right_charge(k);
        const_iterator it = B.basis().left_lower_bound(ar);
        for ( ; it != B_end && it->lc == ar; ++it)
        {
            if (!A_right_basis.has(it->rc)) 
                continue; // trim
            if (!ret.has(A_basis.left_charge(k), it->rc))
                ret.insert(typename DualIndex<SymmGroup>::value_type(A_basis.left_charge(k), it->rc, A_basis.left_size(k), it->rs));
        }
    }
    return ret;
}

/**
 * @brief Returns the index associated with a given b value of the BoundaryMPSProduct object.
 * 
 * The first step of the contraction of a boundary with the MPS sitting on the next site (the 
 * same holds true also for the site_hamiltonian) is the creation of the BoundaryMPSProduct
 * object, that stores the partial contraction of the boundary with the MPS.
 * In some cases, to avoid storing large objects, QCMaquis postpones this contraction 
 * to a later step, specifically to the [Kernel] call.
 * The Kernel call contracts the BoundaryMPSProduct with the MPO and, as a first step,
 * allocates the memory that will accomodate the result of the Kernel().
 * This method is used to estimate the memory required for this allocation step.
 * 
 * @param boundary Input (left) boundary object.
 * @param mult_mps BoundaryMPSProduct object.
 * @param mpo Input Matrix Product Operator
 * @param mps_basis basis of the original mps
 * @param refBasis reference basis for (left) trimming.
 * @param b index of the BoundaryMPSProduct object
 * @param isHermitian if true, activates the Hermitian treatment.
 * @return DualIndex<SymmGroup> Index of the output block_matrix.
 */
template<class Matrix, class OtherMatrix, class SymmGroup>
DualIndex<SymmGroup> T_basis_left(Boundary<OtherMatrix, SymmGroup> const & boundary,
                                  common::BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemms> const & mult_mps,
                                  MPOTensor<Matrix, SymmGroup> const & mpo, DualIndex<SymmGroup> const & mps_basis,
                                  DualIndex<SymmGroup> const & refBasis, typename MPOTensor<Matrix, SymmGroup>::index_type b,
                                  bool isHermitian)
{
    if (mpo.num_row_non_zeros(b) == 1) {
        if (mpo.herm_info.left_skip(b) && isHermitian) {
            return gemm_trim_left_basis(conjugate(boundary[mpo.herm_info.left_conj(b)]), mps_basis, refBasis);
        }
        else {
            return gemm_trim_left_basis(transpose(boundary[b]), mps_basis, refBasis);
        }
    }
    else {
        return mult_mps[b].basis();
    }
}

/**
 * @brief Right counterpart of [T_basis_left]
 * @param boundary Input (right) boundary
 * @param mult_mps MPSBoundaryProduct object
 * @param mpo Input MPO
 * @param mps_basis Basis of the MPS
 * @param refBasis basis used as a reference to trim
 * @param b index of the output MPSBoundaryProduct object.
 * @param isHermitian if true, activates the Hermitian treatment.
 * @return DualIndex<SymmGroup> Index storing the dimensions of the output block_matrix
 */
template<class Matrix, class OtherMatrix, class SymmGroup>
DualIndex<SymmGroup> T_basis_right(Boundary<OtherMatrix, SymmGroup> const & boundary,
                                   common::MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, Gemms> const & mult_mps,
                                   MPOTensor<Matrix, SymmGroup> const & mpo,
                                   DualIndex<SymmGroup> const & mps_basis, DualIndex<SymmGroup> const & refBasis,
                                   typename MPOTensor<Matrix, SymmGroup>::index_type b, bool isHermitian)
{
    if (mpo.num_col_non_zeros(b) == 1)
        if (mpo.herm_info.right_skip(b) && isHermitian)
            return gemm_trim_right_basis(mps_basis, transpose(conjugate(boundary[mpo.herm_info.right_conj(b)])), refBasis);
        else
            return gemm_trim_right_basis(mps_basis, boundary[b], refBasis);
    else
        return mult_mps[b].basis();
}

}
} // namespace abelian
} // namespace contraction

#endif
