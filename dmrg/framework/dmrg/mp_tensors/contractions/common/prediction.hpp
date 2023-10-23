/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef ENGINE_COMMON_PREDICTION_H
#define ENGINE_COMMON_PREDICTION_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/indexing.h"

namespace contraction {
namespace common {

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
static std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
predict_new_state_l2r_sweep(MPSTensor<Matrix, SymmGroup> const & mps, MPOTensor<Matrix, SymmGroup> const & mpo,
                            Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                            double alpha, double cutoff, std::size_t Mmax, bool doPerturbDM, bool verbose=false)
{
    mps.make_left_paired();
    block_matrix<Matrix, SymmGroup> U, V, dm;
    block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
    truncation_results trunc;
    MPSTensor<Matrix, SymmGroup> ret = mps;
    if (doPerturbDM) {
        typename Gemm::gemm()(mps.data(), transpose(conjugate(mps.data())), dm);
        assert( weak_equal(dm.left_basis(), mps.data().left_basis()) );
        Boundary<Matrix, SymmGroup> half_dm
          = left_boundary_tensor_mpo<Matrix, OtherMatrix, SymmGroup, Gemm, Kernel>(mps, left, mpo);
        for (int b = 0; b < half_dm.aux_dim(); ++b) {
            block_matrix<Matrix, SymmGroup> tdm;
            typename Gemm::gemm()(half_dm[b], transpose(conjugate(half_dm[b])), tdm);
            tdm *= alpha;
            for (int k = 0; k < tdm.n_blocks(); ++k) {
                if (mps.data().basis().has(tdm.basis().left_charge(k), tdm.basis().right_charge(k)))
                    dm.match_and_add_block(tdm[k],
                                           tdm.basis().left_charge(k),
                                           tdm.basis().right_charge(k));
            }
        }
        trunc = heev_truncate(dm, U, S, cutoff, Mmax, verbose);
    }
    else {
        block_matrix<Matrix, SymmGroup> tmp = mps.data();
        trunc = svd_truncate(tmp, U, V, S, cutoff, Mmax, verbose);
        //typename Gemm::gemm()(mps.data(), transpose(conjugate(mps.data())), dm);
        //trunc = heev_truncate(dm, U, S, cutoff, Mmax);
    }
    ret.replace_left_paired(U);
    return std::make_pair(ret, trunc);
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
static block_matrix<Matrix, SymmGroup>
getZeroSiteTensorL2R(MPSTensor<Matrix, SymmGroup> mpsNextSite, const MPSTensor<Matrix, SymmGroup>& mpsCurrentSite,
                     const MPSTensor<Matrix, SymmGroup>& U)
{
    mpsCurrentSite.make_left_paired();
    U.make_left_paired();
    block_matrix<Matrix, SymmGroup> tmp;
    // The factor to be included in the following site is SU^T = U^TM
    typename Gemm::gemm()(transpose(conjugate(U.data())), mpsCurrentSite.data(), tmp);
    return tmp;
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
static MPSTensor<Matrix, SymmGroup>
predict_lanczos_l2r_sweep(MPSTensor<Matrix, SymmGroup> mpsNextSite,
                          const MPSTensor<Matrix, SymmGroup>& mpsCurrentSite,
                          const MPSTensor<Matrix, SymmGroup>& U)
{
    auto tmp = getZeroSiteTensorL2R(mpsNextSite, mpsCurrentSite, U);
    mpsNextSite.multiply_from_left(tmp);
    return mpsNextSite;
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
static std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
predict_new_state_r2l_sweep(MPSTensor<Matrix, SymmGroup> const & mps, MPOTensor<Matrix, SymmGroup> const & mpo,
                            Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                            double alpha, double cutoff, std::size_t Mmax, bool doPerturbDM, bool verbose=false)
{
    // Initialization
    mps.make_right_paired();
    block_matrix<Matrix, SymmGroup> dm;
    MPSTensor<Matrix, SymmGroup> ret = mps;
    truncation_results trunc;
    block_matrix<Matrix, SymmGroup> U, V;
    block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
    // Compute "real" density matrix
    mps.make_right_paired();
    if (doPerturbDM) {
        // Add noise
        typename Gemm::gemm()(transpose(conjugate(mps.data())), mps.data(), dm);
        Boundary<Matrix, SymmGroup> half_dm
          = right_boundary_tensor_mpo<Matrix, OtherMatrix, SymmGroup, Gemm, Kernel>(mps, right, mpo);
        for (std::size_t b = 0; b < half_dm.aux_dim(); ++b) {
            block_matrix<Matrix, SymmGroup> tdm;
            typename Gemm::gemm()(transpose(conjugate(half_dm[b])), half_dm[b], tdm);
            tdm *= alpha;
            for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                if (mps.data().basis().has(tdm.basis().left_charge(k), tdm.basis().right_charge(k)))
                    dm.match_and_add_block(tdm[k],
                                           tdm.basis().left_charge(k),
                                           tdm.basis().right_charge(k));
            }
        }
        mps.make_right_paired();
        assert( weak_equal(dm.right_basis(), mps.data().right_basis()) );
        trunc = heev_truncate(dm, U, S, cutoff, Mmax, verbose);
    }
    else {
        block_matrix<Matrix, SymmGroup> tmp = transpose(conjugate(mps.data()));
        trunc = svd_truncate(tmp, U, V, S, cutoff, Mmax, verbose);
    }
    ret.replace_right_paired(adjoint(U));
    return std::make_pair(ret, trunc);
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
static block_matrix<Matrix, SymmGroup>
getZeroSiteTensorR2L(MPSTensor<Matrix, SymmGroup> mpsPreviousSite,
                     MPSTensor<Matrix, SymmGroup> const & mpsCurrentSite,
                     MPSTensor<Matrix, SymmGroup> const & U)
{
    mpsCurrentSite.make_right_paired();
    U.make_right_paired();
    block_matrix<Matrix, SymmGroup> tmp;
    typename Gemm::gemm()(mpsCurrentSite.data(), transpose(conjugate(U.data())), tmp);
    return tmp;
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
static MPSTensor<Matrix, SymmGroup>
predict_lanczos_r2l_sweep(MPSTensor<Matrix, SymmGroup> mpsPreviousSite,
                          const MPSTensor<Matrix, SymmGroup>& mpsCurrentSite,
                          MPSTensor<Matrix, SymmGroup> const & U)
{
    auto tmp = getZeroSiteTensorR2L(mpsPreviousSite, mpsCurrentSite, U);
    mpsPreviousSite.multiply_from_right(tmp);
    return mpsPreviousSite;
}

} // namespace common
} // namespace contraction

#endif
