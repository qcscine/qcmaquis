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

#ifndef SU2_ENGINE_HPP
#define SU2_ENGINE_HPP

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

#include "dmrg/mp_tensors/contractions/common/common.h"

#include "dmrg/mp_tensors/contractions/non-abelian/apply_op.hpp"
#include "dmrg/mp_tensors/contractions/non-abelian/apply_op_rp.hpp"
#include "dmrg/mp_tensors/contractions/non-abelian/gemm.hpp"
#include "dmrg/mp_tensors/contractions/non-abelian/functors.h"
#include "dmrg/mp_tensors/contractions/non-abelian/h_diag.hpp"

namespace contraction {

    using ::contraction::common::BoundaryMPSProduct;
    using ::contraction::common::MPSBoundaryProduct;
    using ::contraction::common::MPSBoundaryProductIndices;

    template <class Matrix, class OtherMatrix, class SymmGroup>
    class Engine<Matrix, OtherMatrix, SymmGroup, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type>
    {
        struct lbtm_functor
        {
            void operator()(size_t b2,
                            contraction::ContractionGrid<Matrix, SymmGroup>& contr_grid,
                            Boundary<OtherMatrix, SymmGroup> const & left,
                            BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> const & left_mult_mps,
                            MPOTensor<Matrix, SymmGroup> const & mpo,
                            DualIndex<SymmGroup> const & ket_basis,
                            Index<SymmGroup> const & right_i,
                            Index<SymmGroup> const & out_left_i,
                            ProductBasis<SymmGroup> const & in_right_pb,
                            ProductBasis<SymmGroup> const & out_left_pb)
            {
                return SU2::lbtm_kernel(b2, contr_grid, left, left_mult_mps, mpo, ket_basis, right_i, out_left_i, in_right_pb, out_left_pb);
            }
        };

        struct rbtm_functor
        {
            void operator()(size_t b1,
                            block_matrix<Matrix, SymmGroup> & ret,
                            Boundary<OtherMatrix, SymmGroup> const & right,
                            MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> const & right_mult_mps,
                            MPOTensor<Matrix, SymmGroup> const & mpo,
                            DualIndex<SymmGroup> const & ket_basis,
                            Index<SymmGroup> const & left_i,
                            Index<SymmGroup> const & out_right_i,
                            ProductBasis<SymmGroup> const & in_left_pb,
                            ProductBasis<SymmGroup> const & out_right_pb)
            {
                return SU2::rbtm_kernel(b1, ret, right, right_mult_mps, mpo, ket_basis, left_i, out_right_i, in_left_pb, out_right_pb);
            }
        };

    public:

        typedef typename contraction::common::Schedule<Matrix, SymmGroup>::schedule_t schedule_t;

        static block_matrix<OtherMatrix, SymmGroup>
        overlap_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          block_matrix<OtherMatrix, SymmGroup> const & left,
                          block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
        {
            return common::overlap_left_step<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms>
                   (bra_tensor, ket_tensor, left, localop);
        }

        static block_matrix<OtherMatrix, SymmGroup>
        overlap_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                           MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                           block_matrix<OtherMatrix, SymmGroup> const & right,
                           block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
        {
            return common::overlap_right_step<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms>
                   (bra_tensor, ket_tensor, right, localop);
        }

        static Boundary<Matrix, SymmGroup>
        left_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                                 Boundary<OtherMatrix, SymmGroup> const & left,
                                 MPOTensor<Matrix, SymmGroup> const & mpo,
                                 Index<SymmGroup> const * in_low = NULL)
        {
            return common::left_boundary_tensor_mpo<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms, lbtm_functor>
                   (mps, left, mpo, in_low);
        }

        static Boundary<Matrix, SymmGroup>
        right_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                                  Boundary<OtherMatrix, SymmGroup> const & right,
                                  MPOTensor<Matrix, SymmGroup> const & mpo,
                                  Index<SymmGroup> const * in_low = NULL)
        {
            return common::right_boundary_tensor_mpo<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms, rbtm_functor>
                   (mps, right, mpo, in_low);
        }

        static Boundary<OtherMatrix, SymmGroup>
        overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                              MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                              Boundary<OtherMatrix, SymmGroup> const & left,
                              MPOTensor<Matrix, SymmGroup> const & mpo)
        {
            return common::overlap_mpo_left_step<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms, lbtm_functor>
                   (bra_tensor, ket_tensor, left, mpo);
        }

        static Boundary<OtherMatrix, SymmGroup>
        overlap_mpo_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                               MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                               Boundary<OtherMatrix, SymmGroup> const & right,
                               MPOTensor<Matrix, SymmGroup> const & mpo)
        {
            return common::overlap_mpo_right_step<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms, rbtm_functor>
                   (bra_tensor, ket_tensor, right, mpo);
        }

        static schedule_t
        right_contraction_schedule(MPSTensor<Matrix, SymmGroup> const & mps,
                                   Boundary<OtherMatrix, SymmGroup> const & right,
                                   MPOTensor<Matrix, SymmGroup> const & mpo)
        {
            return common::create_contraction_schedule(mps, right, mpo, SU2::rbtm_tasks<Matrix, OtherMatrix, SymmGroup>);
        }

        // Single-site prediction
        static std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
        predict_new_state_l2r_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                    MPOTensor<Matrix, SymmGroup> const & mpo,
                                    Boundary<OtherMatrix, SymmGroup> const & left,
                                    Boundary<OtherMatrix, SymmGroup> const & right,
                                    double alpha,
                                    double cutoff,
                                    std::size_t Mmax,
                                    std::size_t Mval
                                   )
        {
            return common::predict_new_state_l2r_sweep<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms, lbtm_functor>
                   (mps, mpo, left, right, alpha, cutoff, Mmax, Mval);
        }

        static MPSTensor<Matrix, SymmGroup>
        predict_lanczos_l2r_sweep(MPSTensor<Matrix, SymmGroup> B,
                                  MPSTensor<Matrix, SymmGroup> const & psi,
                                  MPSTensor<Matrix, SymmGroup> const & A)
        {
            return common::predict_lanczos_l2r_sweep<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms>(B, psi, A);
        }

        static std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
        predict_new_state_r2l_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                    MPOTensor<Matrix, SymmGroup> const & mpo,
                                    Boundary<OtherMatrix, SymmGroup> const & left,
                                    Boundary<OtherMatrix, SymmGroup> const & right,
                                    double alpha,
                                    double cutoff,
                                    std::size_t Mmax,
                                    std::size_t Mval)
        {
            return common::predict_new_state_r2l_sweep<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms, rbtm_functor>
                   (mps, mpo, left, right, alpha, cutoff, Mmax, Mval);
        }

        static MPSTensor<Matrix, SymmGroup>
        predict_lanczos_r2l_sweep(MPSTensor<Matrix, SymmGroup> B,
                                  MPSTensor<Matrix, SymmGroup> const & psi,
                                  MPSTensor<Matrix, SymmGroup> const & A)
        {
            return common::predict_lanczos_r2l_sweep<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms>(B, psi, A);
        }

        static MPSTensor<Matrix, SymmGroup>
        site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
                    Boundary<OtherMatrix, SymmGroup> const & left,
                    Boundary<OtherMatrix, SymmGroup> const & right,
                    MPOTensor<Matrix, SymmGroup> const & mpo,
                    std::vector<SU2::task_capsule<Matrix, SymmGroup> > const & tasks);
    };

} // namespace contraction

#include "dmrg/mp_tensors/contractions/non-abelian/site_hamil.hpp"

#endif
