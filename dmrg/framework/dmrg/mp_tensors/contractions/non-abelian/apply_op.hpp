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

#ifndef CONTRACTIONS_SU2_APPLY_OP_HPP
#define CONTRACTIONS_SU2_APPLY_OP_HPP

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/contractions/non-abelian/gsl_coupling.h"

namespace SU2 {

    inline
    double mod_coupling(int two_ja, int two_jb, int two_jc,
                        int two_jd, int two_je, int two_jf,
                        int two_jg, int two_jh, int two_ji)
    {
        return sqrt( (two_jg+1.) * (two_jh+1.) * (two_jc+1.) * (two_jf+1.) ) *
               gsl_sf_coupling_9j(two_ja, two_jb, two_jc,
                                  two_jd, two_je, two_jf,
                                  two_jg, two_jh, two_ji);
    }
}

namespace contraction {
namespace SU2 {

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel(size_t b2,
                     ContractionGrid<Matrix, SymmGroup>& contr_grid,
                     Boundary<OtherMatrix, SymmGroup> const & left,
                     std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     DualIndex<SymmGroup> const & ket_basis,
                     Index<SymmGroup> const & right_i,
                     Index<SymmGroup> const & out_left_i,
                     ProductBasis<SymmGroup> const & in_right_pb,
                     ProductBasis<SymmGroup> const & out_left_pb)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename SymmGroup::charge charge;

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();

            block_matrix<Matrix, SymmGroup> const & T = left_mult_mps[b1];
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            block_matrix<Matrix, SymmGroup> const & W = access.op;
            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);

            ret.spin = couple(left[b1].spin, W.spin);
            int a = left[b1].spin.get(), k = W.spin.get(), ap = ret.spin.get();

            if (k != W.spin.get()) { maquis::cout << "k " << k << "W.spin.get() " << W.spin.get() << std::endl; }

            for (size_t lblock = 0; lblock < left[b1].n_blocks(); ++lblock) {

                charge lc = left[b1].basis().right_charge(lblock); // left comes as left^T !
                charge mc = left[b1].basis().left_charge(lblock);

                std::pair<const_iterator, const_iterator>
                  er = std::equal_range(ket_basis.begin(), ket_basis.end(),
                    dual_index_detail::QnBlock<SymmGroup>(mc, SymmGroup::IdentityCharge, 0, 0), dual_index_detail::gt_row<SymmGroup>());

                for (const_iterator it = er.first; it != er.second; ++it)
                {
                    size_t matched_block = std::distance(ket_basis.begin(), it);
                    charge rc = ket_basis[matched_block].rc;
                    size_t t_block = T.basis().position(lc, rc); // t_block != lblock in general

                    for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                    {
                        charge phys_in = W.basis().left_charge(w_block);
                        charge phys_out = W.basis().right_charge(w_block);

                        charge out_r_charge = SymmGroup::fuse(rc, phys_in);
                        if (!right_i.has(out_r_charge)) continue;

                        charge out_l_charge = SymmGroup::fuse(lc, phys_out);
                        if (!right_i.has(out_l_charge) || !out_left_i.has(out_l_charge)) continue;

                        size_t r_size = right_i.size_of_block(out_r_charge);

                        size_t o = ret.find_block(out_l_charge, out_r_charge);
                        if ( o == ret.n_blocks() ) {
                            o = ret.insert_block(Matrix(1,1), out_l_charge, out_r_charge);
                            ret.resize_block(o, out_left_i.size_of_block(out_l_charge), r_size);
                        }

                        int i  = lc[1], ip = out_l_charge[1];
                        int j  = mc[1], jp  = out_r_charge[1];
                        int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                        //typename Matrix::value_type coupling_coeff = ::SU2::mod_coupling(j, two_s, jp, left[b1].spin.get(),W.spin.get(),ret.spin.get(), i, two_sp, ip);
                        typename Matrix::value_type coupling_coeff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                        if (std::abs(coupling_coeff) < 1.e-40) continue;
                        coupling_coeff *= sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.))) * access.scale;

                        size_t phys_s1 = W.basis().left_size(w_block);
                        size_t phys_s2 = W.basis().right_size(w_block);
                        size_t in_right_offset = in_right_pb(phys_in, out_r_charge);
                        size_t out_left_offset = out_left_pb(phys_out, lc);
                        Matrix const & wblock = W[w_block];
                        Matrix const & iblock = T[t_block];
                        Matrix & oblock = ret[o];

                        //maquis::cout << "access " << mc << " + " << phys_in<< "|" << out_r_charge << " -- "
                        //         << lc << " + " << phys_out << "|" << out_l_charge
                        //         << " T(" << lc << "," << rc << "): +" << in_right_offset
                        //         << "(" << T.basis().left_size(t_block) << "x" << r_size << ")|" << T.basis().right_size(t_block) << " -> +"
                        //         << out_left_offset << "(" << T.basis().left_size(t_block) << "x" << r_size << ")|" << out_left_i.size_of_block(out_l_charge)
                        //         << " * " << j << two_s << jp << a << k << ap << i << two_sp << ip << " " << coupling_coeff * wblock(0,0) << std::endl;

                        maquis::dmrg::detail::lb_tensor_mpo(oblock, iblock, wblock,
                                out_left_offset, in_right_offset,
                                phys_s1, phys_s2, T.basis().left_size(t_block), r_size, coupling_coeff);
                    }
                }
            }
        } // b1
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void rbtm_kernel(size_t b1,
                block_matrix<Matrix, SymmGroup> & ret,
                Boundary<OtherMatrix, SymmGroup> const & right,
                std::vector<block_matrix<Matrix, SymmGroup> > const & right_mult_mps,
                MPOTensor<Matrix, SymmGroup> const & mpo,
                DualIndex<SymmGroup> const & ket_basis,
                Index<SymmGroup> const & left_i,
                Index<SymmGroup> const & out_right_i,
                ProductBasis<SymmGroup> const & in_left_pb,
                ProductBasis<SymmGroup> const & out_right_pb)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename SymmGroup::charge charge;

        row_proxy row_b1 = mpo.row(b1);
        for (typename row_proxy::const_iterator row_it = row_b1.begin(); row_it != row_b1.end(); ++row_it) {
            index_type b2 = row_it.index();

            block_matrix<Matrix, SymmGroup> const & T = right_mult_mps[b2];
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            block_matrix<Matrix, SymmGroup> const & W = access.op;

            // note minus sign, as we move the boundary to the left
            ret.spin = couple(right[b2].spin, -W.spin);
            int ap = right[b2].spin.get(), k = W.spin.get(), a = ret.spin.get();

            for (size_t ketblock = 0; ketblock < ket_basis.size(); ++ketblock) {

                charge lc = ket_basis[ketblock].lc;
                charge mc = ket_basis[ketblock].rc;

                std::pair<const_iterator, const_iterator>
                  er = std::equal_range(right[b2].basis().begin(), right[b2].basis().end(),
                    dual_index_detail::QnBlock<SymmGroup>(mc, SymmGroup::IdentityCharge, 0, 0), dual_index_detail::gt_row<SymmGroup>());

                for (const_iterator it = er.first; it != er.second; ++it)
                {
                    size_t matched_block = std::distance(right[b2].basis().begin(), it);
                    charge rc = right[b2].basis()[matched_block].rc;
                    size_t t_block = T.basis().position(lc, rc); // t_block != ketblock in general

                    for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                    {
                        charge phys_in = W.basis().left_charge(w_block);
                        charge phys_out = W.basis().right_charge(w_block);

                        charge out_l_charge = SymmGroup::fuse(lc, -phys_in);
                        if (!left_i.has(out_l_charge)) continue;

                        charge out_r_charge = SymmGroup::fuse(rc, -phys_out);
                        if (!left_i.has(out_r_charge) || !out_right_i.has(out_r_charge)) continue;

                        size_t l_size = left_i.size_of_block(out_l_charge);

                        size_t o = ret.find_block(out_l_charge, out_r_charge);
                        if ( o == ret.n_blocks() ) {
                            o = ret.insert_block(Matrix(1,1), out_l_charge, out_r_charge);
                            ret.resize_block(o, l_size, out_right_i.size_of_block(out_r_charge));
                        }

                        int i = out_r_charge[1], ip = rc[1];
                        int j = out_l_charge[1], jp  = mc[1];
                        int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                        typename Matrix::value_type coupling_coeff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                        if (std::abs(coupling_coeff) < 1.e-40) continue;
                        coupling_coeff *= sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.))) * access.scale;

                        size_t phys_s1 = W.basis().left_size(w_block);
                        size_t phys_s2 = W.basis().right_size(w_block);
                        size_t in_left_offset = in_left_pb(phys_in, out_l_charge);
                        size_t out_right_offset = out_right_pb(phys_out, rc);
                        Matrix const & wblock = W[w_block];
                        Matrix const & iblock = T[t_block];
                        Matrix & oblock = ret[o];

                        //maquis::cout << "access " << mc << " - " << phys_out<< "|" << out_r_charge << " -- "
                        //         << lc << " - " << phys_in << "|" << out_l_charge
                        //         << " T(" << lc << "," << rc << "): +" << in_left_offset
                        //         << "(" << l_size << "x" << T.basis().right_size(t_block) << ")|" << T.basis().left_size(t_block) << " -> +"
                        //         << out_right_offset << "(" << l_size << "x" << T.basis().right_size(t_block) << ")|" << out_right_i.size_of_block(out_r_charge)
                        //         << " * " << j << two_s << jp << a << k << ap << i << two_sp << ip << " " << coupling_coeff * wblock(0,0) << std::endl;

                        maquis::dmrg::detail::rb_tensor_mpo(oblock, iblock, wblock,
                                out_right_offset, in_left_offset,
                                phys_s1, phys_s2, l_size, T.basis().right_size(t_block), coupling_coeff);
                    }
                }
            }
        } // b1
    }
} // namespace SU2
} // namespace contraction

#endif
