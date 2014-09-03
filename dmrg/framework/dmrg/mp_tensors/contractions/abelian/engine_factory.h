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

#ifndef ABELIAN_ENGINE_FACTORY_H
#define ABELIAN_ENGINE_FACTORY_H

#include <boost/shared_ptr.hpp>

#include "dmrg/mp_tensors/contractions/abelian/engine.hpp"

namespace contraction {

    // Forward declaration
    template <class Matrix, class OtherMatrix, class SymmGroup>
    class EngineFactory;
}

namespace contraction {

    template <class Matrix, class OtherMatrix, class SymmGroup>
    class AbelianEngineFactory : public EngineFactory<Matrix, OtherMatrix, SymmGroup>
    {
        typedef boost::shared_ptr<contraction::Engine<Matrix, OtherMatrix, SymmGroup> > engine_ptr;

    public:
        virtual engine_ptr makeEngine() { return engine_ptr(new AbelianEngine<Matrix, OtherMatrix, SymmGroup>()); }

        static block_matrix<OtherMatrix, SymmGroup>
        overlap_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          block_matrix<OtherMatrix, SymmGroup> const & left,
                          block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template overlap_left_step<AbelianGemms>(bra_tensor, ket_tensor, left, localop);
        }

        static block_matrix<OtherMatrix, SymmGroup>
        overlap_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                           MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                           block_matrix<OtherMatrix, SymmGroup> const & right,
                           block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template overlap_right_step<AbelianGemms>(bra_tensor, ket_tensor, right, localop);
        }

        static Boundary<Matrix, SymmGroup>
        left_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                                 Boundary<OtherMatrix, SymmGroup> const & left,
                                 MPOTensor<Matrix, SymmGroup> const & mpo,
                                 Index<SymmGroup> const * in_low = NULL)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template left_boundary_tensor_mpo<AbelianGemms, lbtm_functor>
                   (mps, left, mpo, in_low);
        }

        static Boundary<Matrix, SymmGroup>
        right_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                                  Boundary<OtherMatrix, SymmGroup> const & right,
                                  MPOTensor<Matrix, SymmGroup> const & mpo,
                                  Index<SymmGroup> const * in_low = NULL)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template right_boundary_tensor_mpo<AbelianGemms, rbtm_functor>
                   (mps, right, mpo, in_low);
        }

        static Boundary<OtherMatrix, SymmGroup>
        overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                              MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                              Boundary<OtherMatrix, SymmGroup> const & left,
                              MPOTensor<Matrix, SymmGroup> const & mpo)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template overlap_mpo_left_step<AbelianGemms, lbtm_functor>
                   (bra_tensor, ket_tensor, left, mpo);
        }

        static Boundary<OtherMatrix, SymmGroup>
        overlap_mpo_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                               MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                               Boundary<OtherMatrix, SymmGroup> const & right,
                               MPOTensor<Matrix, SymmGroup> const & mpo)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template overlap_mpo_right_step<AbelianGemms, rbtm_functor>
                   (bra_tensor, ket_tensor, right, mpo);
        }

        static MPSTensor<Matrix, SymmGroup>
        site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
                    Boundary<OtherMatrix, SymmGroup> const & left,
                    Boundary<OtherMatrix, SymmGroup> const & right,
                    MPOTensor<Matrix, SymmGroup> const & mpo);

        static std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
        predict_new_state_l2r_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                    MPOTensor<Matrix, SymmGroup> const & mpo,
                                    Boundary<OtherMatrix, SymmGroup> const & left,
                                    Boundary<OtherMatrix, SymmGroup> const & right,
                                    double alpha, double cutoff, std::size_t Mmax)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template predict_new_state_l2r_sweep<AbelianGemms, lbtm_functor>
                  (mps, mpo, left, right, alpha, cutoff, Mmax);
        }

        static MPSTensor<Matrix, SymmGroup>
        predict_lanczos_l2r_sweep(MPSTensor<Matrix, SymmGroup> B,
                                  MPSTensor<Matrix, SymmGroup> const & psi,
                                  MPSTensor<Matrix, SymmGroup> const & A)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template predict_lanczos_l2r_sweep<AbelianGemms>(B, psi, A);
        }

        static std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
        predict_new_state_r2l_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                    MPOTensor<Matrix, SymmGroup> const & mpo,
                                    Boundary<OtherMatrix, SymmGroup> const & left,
                                    Boundary<OtherMatrix, SymmGroup> const & right,
                                    double alpha, double cutoff, std::size_t Mmax)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template predict_new_state_r2l_sweep<AbelianGemms, rbtm_functor>
                   (mps, mpo, left, right, alpha, cutoff, Mmax);
        }

        static MPSTensor<Matrix, SymmGroup>
        predict_lanczos_r2l_sweep(MPSTensor<Matrix, SymmGroup> B,
                                  MPSTensor<Matrix, SymmGroup> const & psi,
                                  MPSTensor<Matrix, SymmGroup> const & A)
        {
            return EngineFactory<Matrix, OtherMatrix, SymmGroup>::template predict_lanczos_r2l_sweep<AbelianGemms>(B, psi, A);
        }

    };

} // namespace contraction

#endif
