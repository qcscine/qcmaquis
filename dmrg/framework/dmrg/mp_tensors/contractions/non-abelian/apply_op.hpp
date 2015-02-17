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

#include "dmrg/block_matrix/symmetry/gsl_coupling.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

namespace contraction {
namespace SU2 {

    template<class Matrix, class SymmGroup>
    bool track_it(block_matrix<Matrix, SymmGroup> const & op)
    {
        if (op.n_blocks() == 1) return true;
        if (std::abs(op.basis().left_charge(0)[0] - op.basis().right_charge(0)[0]) == 2)
            return true;
        return false;
    }

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

        bool debug = false;
        //bool debug = (mpo.col_dim() == 1);
        //bool debug = track_it(W);
        if(debug) {
            maquis::cout << "right_i   : " << right_i << std::endl;
            maquis::cout << "out_left_i: " << out_left_i << std::endl;
            maquis::cout << "ket_basis : " << ket_basis << std::endl << std::endl;
        }

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();

            block_matrix<Matrix, SymmGroup> const & T = left_mult_mps[b1];
            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1,b2);

            if (debug) maquis::cout << "lbtm b1,b2: " << b1 << "," << b2 << std::endl << std::endl;

        for (std::size_t op_index = 0; op_index < access.size(); ++op_index)
        {
            block_matrix<Matrix, SymmGroup> const & W = access.op(op_index);
            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);

            ret.spin = couple(left[b1].spin, W.spin);
            int a = left[b1].spin.get(), k = W.spin.get(), ap = ret.spin.get();
            assert(a == mpo.left_spin(b1).get() && ap == mpo.right_spin(b2).get());
            if (debug) {maquis::cout << "  WBase: " << W.basis() << std::endl;}

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

                        int i = lc[1], ip = out_l_charge[1];
                        int j = mc[1], jp = out_r_charge[1];
                        int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                        //typename Matrix::value_type coupling_coeff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip, debug);
                        //if (std::abs(coupling_coeff) < 1.e-40) continue;
                        //coupling_coeff *= sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.))) * access.scale;

                        size_t phys_s1 = W.basis().left_size(w_block);
                        size_t phys_s2 = W.basis().right_size(w_block);
                        size_t in_right_offset = in_right_pb(phys_in, out_r_charge);
                        size_t out_left_offset = out_left_pb(phys_out, lc);
                        Matrix const & wblock = W[w_block];
                        Matrix const & iblock = T[t_block];
                        Matrix & oblock = ret[o];

                        if(debug)
                        maquis::cout << "    access " << mc << " + " << phys_in<< "|" << out_r_charge << " -- "
                                 << lc << " + " << phys_out << "|" << out_l_charge
                                 << " T(" << lc << "," << rc << "): +" << in_right_offset
                                 << "(" << T.basis().left_size(t_block) << "x" << r_size << ")|" << T.basis().right_size(t_block) << " -> +"
                                 << out_left_offset << "(" << T.basis().left_size(t_block) << "x" << r_size << ")|" << out_left_i.size_of_block(out_l_charge)
                                 << std::endl;
                                 //<< " * " << j << two_s << jp << a << k << ap << i << two_sp << ip << " " << coupling_coeff * wblock(0,0) << std::endl;

                        //maquis::dmrg::detail::lb_tensor_mpo(oblock, iblock, wblock,
                        //        out_left_offset, in_right_offset,
                        //        phys_s1, phys_s2, T.basis().left_size(t_block), r_size, coupling_coeff);

                        typename Matrix::value_type c_eff = 555555;//coupling_coeff;
                        size_t ldim = T.basis().left_size(t_block);
                        for(size_t rr = 0; rr < r_size; ++rr) {
                            for(size_t ss1 = 0; ss1 < phys_s1; ++ss1) {
                                for(size_t ss2 = 0; ss2 < phys_s2; ++ss2) {
                                    if (ss1 == 2 && ss2 == 2) {
                                        int ip_loc = (ip == 0 && k == 2) ? ip : ip;
                                        int jp_loc = (jp == 0 && k == 2) ? jp : jp;
                                        if(debug) maquis::cout << "      " << std::setw(12) << wblock(ss1,ss2);
                                        c_eff = ::SU2::mod_coupling(j, 2, jp_loc, a,k,ap, i, 2, ip_loc, debug);
                                    }
                                    else if (ss1 == 2) {
                                        // jp ~ out_r_charge ~ phys_in
                                        int jp_loc = (jp == 0 && k == 2) ? jp : jp;
                                        if(debug) maquis::cout << "      " << std::setw(12) << wblock(ss1,ss2);
                                        c_eff = ::SU2::mod_coupling(j, 2, jp_loc, a,k,ap, i, two_sp, ip, debug);
                                    }
                                    else if (ss2 == 2) {
                                        // ip ~ out_l_charge ~ phys_out
                                        int ip_loc = (ip == 0 && k == 2) ? ip : ip;
                                        if(debug) maquis::cout << "      " << std::setw(12) << wblock(ss1,ss2);
                                        c_eff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, 2, ip_loc, debug);
                                    }
                                    else {
                                        if(debug) maquis::cout << "      " << std::setw(12) << wblock(ss1,ss2);
                                        c_eff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip, debug);
                                        //c_eff = coupling_coeff;
                                    }

                                    c_eff *= sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.))) * access.scale(op_index);

                                    typename Matrix::value_type alfa_t = wblock(ss1, ss2) * c_eff;
                                    maquis::dmrg::detail::iterator_axpy(&iblock(0, in_right_offset + ss1*r_size + rr),
                                                                        &iblock(0, in_right_offset + ss1*r_size + rr) + ldim, // bugbug
                                                                        &oblock(out_left_offset + ss2*ldim, rr),
                                                                        alfa_t);
                                }
                            }
                        }
                        if(debug) maquis::cout << "\n    " << ret << std::endl;

                    } // wblock
                } // ket matches
            } // lblock
        } // op_index
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
            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1,b2);

        for (std::size_t op_index = 0; op_index < access.size(); ++op_index)
        {
            block_matrix<Matrix, SymmGroup> const & W = access.op(op_index);

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
                        int j = out_l_charge[1], jp = mc[1];
                        int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                        typename Matrix::value_type coupling_coeff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                        if (std::abs(coupling_coeff) < 1.e-40) continue;
                        coupling_coeff *= sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.))) * access.scale(op_index);

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
        } // op_index
        } // b1
    }
} // namespace SU2
} // namespace contraction

#endif
