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

#ifndef CONTRACTIONS_SU2_APPLY_OP_RP_HPP
#define CONTRACTIONS_SU2_APPLY_OP_RP_HPP

#include "dmrg/block_matrix/symmetry/gsl_coupling.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

namespace contraction {
namespace SU2 {

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel_rp(size_t b2,
                        ContractionGrid<Matrix, SymmGroup>& contr_grid,
                        Boundary<OtherMatrix, SymmGroup> const & left,
                        std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
                        MPOTensor<Matrix, SymmGroup> const & mpo,
                        MPSTensor<Matrix, SymmGroup> const & mps,
                        Index<SymmGroup> const & right_i,
                        Index<SymmGroup> const & out_left_i,
                        ProductBasis<SymmGroup> const & in_right_pb,
                        ProductBasis<SymmGroup> const & out_left_pb,
                        DualIndex<SymmGroup> const & right_basis)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename SymmGroup::charge charge;

        DualIndex<SymmGroup> const & ket_basis = mps.data().basis();

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();

            block_matrix<Matrix, SymmGroup> local_T;
            block_matrix<Matrix, SymmGroup> const * Tp = &local_T;
            if (mpo.num_row_non_zeros(b1) == 1)
                ::SU2::gemm_trim_left(transpose(left[b1]), mps.data(), local_T);
            else
                Tp = &left_mult_mps[b1];

            block_matrix<Matrix, SymmGroup> const & T = *Tp;

            //block_matrix<Matrix, SymmGroup> const & T = left_mult_mps[b1];
            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1,b2);

        for (std::size_t op_index = 0; op_index < access.size(); ++op_index)
        {
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op(op_index);
            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);

            int a = mpo.left_spin(b1).get(), k = W.spin().get(), ap = mpo.right_spin(b2).get();

            for (size_t lblock = 0; lblock < left[b1].n_blocks(); ++lblock) {

                charge lc = left[b1].basis().right_charge(lblock); // left comes as left^T !
                charge mc = left[b1].basis().left_charge(lblock);

                if (!ket_basis.left_has(lc)) continue; // ket_basis needs to be left_paired

                const_iterator it = ket_basis.left_lower_bound(mc);
                for ( ; it != ket_basis.end() && it->lc == mc; ++it)
                {
                    charge rc = it->rc;
                    size_t t_block = T.basis().position(lc, rc); // t_block != lblock in general

                    for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                    {
                        charge phys_in = W.basis().left_charge(w_block);
                        charge phys_out = W.basis().right_charge(w_block);

                        charge out_r_charge = SymmGroup::fuse(rc, phys_in);
                        if (!right_i.has(out_r_charge)) continue;

                        charge out_l_charge = SymmGroup::fuse(lc, phys_out);
                        if (!right_i.has(out_l_charge)) continue; // can also probe out_left_i, but right_i has the same charges

                        // optimization: avoid creating blocks which will not be used
                        const_iterator right_it = right_basis.left_lower_bound(out_r_charge);
                        bool right_match = false;
                        for ( ; right_it != right_basis.end() && right_it->lc == out_r_charge; ++right_it)
                        {
                            if (out_l_charge == right_it->rc)
                                right_match = true;
                        }
                        if (!right_match) continue;
                        //////////////////////////////////////////

                        charge out_r_charge_rp = SymmGroup::fuse(out_r_charge, -phys_out);

                        size_t r_size = right_i.size_of_block(out_r_charge);

                        size_t o = ret.find_block(lc, out_r_charge_rp);
                        if ( o == ret.n_blocks() ) {
                            o = ret.insert_block(Matrix(T.basis().left_size(t_block),
                                                 in_right_pb.size(-phys_out, out_r_charge)), lc, out_r_charge_rp);
                        }

                        int i = SymmGroup::spin(lc), ip = SymmGroup::spin(out_l_charge);
                        int j = SymmGroup::spin(mc), jp = SymmGroup::spin(out_r_charge);
                        int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                        //typename Matrix::value_type coupling_coeff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                        //if (std::abs(coupling_coeff) < 1.e-40) continue;
                        //coupling_coeff *= sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.))) * access.scale(op_index);

                        typename Matrix::value_type prefactor = sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.))) * access.scale(op_index);
                        typename Matrix::value_type couplings[4];
                        couplings[0] = prefactor * ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                        couplings[1] = prefactor * ::SU2::mod_coupling(j, 2,     jp, a,k,ap, i, two_sp, ip);
                        couplings[2] = prefactor * ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, 2,      ip);
                        couplings[3] = prefactor * ::SU2::mod_coupling(j, 2,     jp, a,k,ap, i, 2,      ip);

                        size_t phys_s1 = W.basis().left_size(w_block);
                        size_t phys_s2 = W.basis().right_size(w_block);
                        size_t in_right_offset = in_right_pb(phys_in, out_r_charge);
                        size_t out_right_offset = in_right_pb(phys_out, out_r_charge);
                        Matrix const & wblock = W[w_block];
                        Matrix const & iblock = T[t_block];
                        Matrix & oblock = ret[o];

                        //maquis::dmrg::detail::lb_tensor_mpo(oblock, iblock, wblock,
                        //        out_left_offset, in_right_offset,
                        //        phys_s1, phys_s2, T.basis().left_size(t_block), r_size, coupling_coeff);

                        typedef typename SparseOperator<Matrix, SymmGroup>::const_iterator block_iterator;
                        std::pair<block_iterator, block_iterator> blocks = W.get_sparse().block(w_block);

                        size_t ldim = T.basis().left_size(t_block);

                        const size_t chunk = 1024;
                        const size_t blength = r_size*ldim;
                        for(size_t rr = 0; rr < blength/chunk; ++rr) {
                            for( block_iterator it = blocks.first; it != blocks.second; ++it)
                            {
                                std::size_t ss1 = it->row;
                                std::size_t ss2 = it->col;
                                std::size_t rspin = it->row_spin;
                                std::size_t cspin = it->col_spin;
                                std::size_t casenr = 0;
                                if (rspin == 2 && cspin == 2) casenr = 3;
                                else if (rspin == 2) casenr = 1;
                                else if (cspin == 2) casenr = 2;

                                typename Matrix::value_type alfa_t = it->coefficient * couplings[casenr];

                                assert(rr + chunk <= r_size*ldim);
                                maquis::dmrg::detail::iterator_axpy(&iblock(0, in_right_offset + ss1*r_size) + rr*chunk,
                                                                    &iblock(0, in_right_offset + ss1*r_size) + rr*chunk + chunk,
                                                                    &oblock(0, out_right_offset + ss2*r_size) + rr*chunk,
                                                                    alfa_t);
                            }
                        }

                        for( block_iterator it = blocks.first; it != blocks.second; ++it)
                        {
                            std::size_t ss1 = it->row;
                            std::size_t ss2 = it->col;
                            std::size_t rspin = it->row_spin;
                            std::size_t cspin = it->col_spin;
                            std::size_t casenr = 0;
                            if (rspin == 2 && cspin == 2) casenr = 3;
                            else if (rspin == 2) casenr = 1;
                            else if (cspin == 2) casenr = 2;

                            typename Matrix::value_type alfa_t = it->coefficient * couplings[casenr];

                            std::size_t start = blength - blength%chunk;
                            maquis::dmrg::detail::iterator_axpy(&iblock(0, in_right_offset + ss1*r_size) + start,
                                                                &iblock(0, in_right_offset + ss1*r_size) + blength,
                                                                &oblock(0, out_right_offset + ss2*r_size) + start,
                                                                alfa_t);
                        }

                    } // wblock
                } // ket matches
            } // lblock
        } // op_index
        } // b1
    }

} // namespace SU2
} // namespace contraction

#endif
