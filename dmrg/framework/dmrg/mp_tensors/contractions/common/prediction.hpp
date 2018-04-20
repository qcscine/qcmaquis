/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef ENGINE_COMMON_PREDICTION_H
#define ENGINE_COMMON_PREDICTION_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/indexing.h"

namespace contraction {
    namespace common {
        //
        // PREDICT_NEW_STATE_L2R_SWEEP
        // ---------------------------
        // Move the sweep procedure one step to the right in the optimization and
        // generates the new MPS.
        template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
        static std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
        predict_new_state_l2r_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                    MPOTensor<Matrix, SymmGroup> const & mpo,
                                    Boundary<OtherMatrix, SymmGroup> const & left,
                                    Boundary<OtherMatrix, SymmGroup> const & right,
                                    double alpha,
                                    double cutoff,
                                    std::size_t Mmax,
                                    const std::vector<size_t>& keeps = std::vector<size_t>())
        {
            mps.make_left_paired();
            block_matrix<Matrix, SymmGroup> dm;
            typename Gemm::gemm()(mps.data(), transpose(conjugate(mps.data())), dm);
            
            Boundary<Matrix, SymmGroup> half_dm
                = left_boundary_tensor_mpo<Matrix, OtherMatrix, SymmGroup, Gemm, Kernel>(mps, left, mpo);
            
            mps.make_left_paired();
            for (std::size_t b = 0; b < half_dm.aux_dim(); ++b)
            {
                block_matrix<Matrix, SymmGroup> tdm;
                typename Gemm::gemm()(half_dm[b], transpose(conjugate(half_dm[b])), tdm);
                
                tdm *= alpha;
                for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                    if (mps.data().basis().has(tdm.basis().left_charge(k), tdm.basis().right_charge(k)))
                        dm.match_and_add_block(tdm[k],
                                               tdm.basis().left_charge(k),
                                               tdm.basis().right_charge(k));
                }
            }

            mps.make_left_paired();
            assert( weak_equal(dm.left_basis(), mps.data().left_basis()) );
            
            block_matrix<Matrix, SymmGroup> U;
            block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
            truncation_results trunc = heev_truncate(dm, U, S, cutoff, Mmax, true, keeps);
            // Prepares results
            MPSTensor<Matrix, SymmGroup> ret = mps;
            ret.replace_left_paired(U);
            return std::make_pair(ret, trunc);
        }
        
                
        //
        // PREDICT_NEW_STATE_R2L_SWEEP
        // ---------------------------
        // Move the sweep procedure one step to the left in the optimization and
        // generates the new MPS. The overall structure is taken from a similar method,
        // predict_new_state_l2r_sweep
        template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
        static std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
        predict_new_state_r2l_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                        MPOTensor<Matrix, SymmGroup> const & mpo,
                                        Boundary<OtherMatrix, SymmGroup> const & left,
                                        Boundary<OtherMatrix, SymmGroup> const & right,
                                        double alpha,
                                        double cutoff,
                                        std::size_t Mmax,
                                        const std::vector<size_t>& keeps = std::vector<size_t>())
        {
            // Initialization
            mps.make_right_paired();
            block_matrix<Matrix, SymmGroup> dm;
            // Compute "real" density matrix
            typename Gemm::gemm()(transpose(conjugate(mps.data())), mps.data(), dm);
            // Add noise
            Boundary<Matrix, SymmGroup> half_dm
                = right_boundary_tensor_mpo<Matrix, OtherMatrix, SymmGroup, Gemm, Kernel>(mps, right, mpo);
            mps.make_right_paired();
            for (std::size_t b = 0; b < half_dm.aux_dim(); ++b)
            {
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
            // Truncate based on the SVD results
            mps.make_right_paired();
            assert( weak_equal(dm.right_basis(), mps.data().right_basis()) );
            block_matrix<Matrix, SymmGroup> U;
            block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
            truncation_results trunc = heev_truncate(dm, U, S, cutoff, Mmax, true, keeps);
            // Prepare output results
            MPSTensor<Matrix, SymmGroup> ret = mps;
            ret.replace_right_paired(adjoint(U));
            return std::make_pair(ret, trunc);
        }

        template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
        static MPSTensor<Matrix, SymmGroup>
        predict_lanczos_l2r_sweep(MPSTensor<Matrix, SymmGroup> B,
                                  MPSTensor<Matrix, SymmGroup> const & psi,
                                  MPSTensor<Matrix, SymmGroup> const & A)
        {
            psi.make_left_paired();
            A.make_left_paired();

            block_matrix<Matrix, SymmGroup> tmp;
            typename Gemm::gemm()(transpose(conjugate(A.data())), psi.data(), tmp);
            B.multiply_from_left(tmp);

            return B;
        }

        template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
        static MPSTensor<Matrix, SymmGroup>
        predict_lanczos_r2l_sweep(MPSTensor<Matrix, SymmGroup> B,
                                  MPSTensor<Matrix, SymmGroup> const & psi,
                                  MPSTensor<Matrix, SymmGroup> const & A)
        {
            psi.make_right_paired();
            A.make_right_paired();
            
            block_matrix<Matrix, SymmGroup> tmp;
            typename Gemm::gemm()(psi.data(), transpose(conjugate(A.data())), tmp);
            
            B.multiply_from_right(tmp);
            
            return B;
        }

    } // namespace common
} // namespace contraction

#endif
