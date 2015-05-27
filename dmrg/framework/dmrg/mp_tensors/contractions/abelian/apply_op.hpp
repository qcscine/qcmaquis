/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef CONTRACTIONS_APPLY_OP_H
#define CONTRACTIONS_APPLY_OP_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/block_matrix/indexing.h"

namespace contraction {
    namespace abelian {

    using ::contraction::ContractionGrid;

    template<class Matrix, class SymmGroup>
    void lbtm_kernel_allocate(size_t b2,
                              ContractionGrid<Matrix, SymmGroup>& contr_grid,
                              std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
                              MPOTensor<Matrix, SymmGroup> const & mpo,
                              Index<SymmGroup> const & right_i,
                              Index<SymmGroup> const & out_left_i)
    {
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::col_proxy col_proxy;

        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();
            block_matrix<Matrix, SymmGroup> const & T = left_mult_mps[b1];                    if(T.n_blocks() == 0) continue;
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op;                            if(W.n_blocks() == 0) continue;

            charge operator_delta = SymmGroup::fuse(W.basis().right_charge(0), -W.basis().left_charge(0));
            charge        T_delta = SymmGroup::fuse(T.basis().right_charge(0), -T.basis().left_charge(0));
            charge    total_delta = SymmGroup::fuse(operator_delta, -T_delta);
        
            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);

            for(size_t r = 0; r < right_i.size(); ++r){
                charge out_r_charge = right_i[r].first;
                charge out_l_charge = SymmGroup::fuse(out_r_charge, total_delta);             if(!out_left_i.has(out_l_charge)) continue;
                size_t r_size = right_i[r].second;
                if(ret.find_block(out_l_charge, out_r_charge) == ret.n_blocks())
                    //ret.resize_block(ret.insert_block(Matrix(1,1), out_l_charge, out_r_charge), 
                    //                 out_left_i.size_of_block(out_l_charge), r_size);
                    ret.insert_block(Matrix(out_left_i.size_of_block(out_l_charge), r_size), out_l_charge, out_r_charge);
            }
        }
        contr_grid.index_sizes(b2);
    }

    // SK: New version which generates same output but uses right-paired input.
    //     The charge delta optimization is indepent from the changes needed to
    //     skip the preceding reshapes.
    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel_execute(size_t b2,
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
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::col_proxy col_proxy;

        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();
            block_matrix<Matrix, SymmGroup> const & T = left_mult_mps[b1];                    if(T.n_blocks() == 0) continue;
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op;                            if(W.n_blocks() == 0) continue;

            // charge deltas are constant for all blocks
            charge operator_delta = SymmGroup::fuse(W.basis().right_charge(0), -W.basis().left_charge(0));
            charge        T_delta = SymmGroup::fuse(T.basis().right_charge(0), -T.basis().left_charge(0));
            charge    total_delta = SymmGroup::fuse(operator_delta, -T_delta);
        
            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);
            parallel::guard group(contr_grid.where(b1,b2), contr_grid.granularity);
            parallel::scheduler_size_indexed scheduler(ret);

            for (size_t r = 0; r < right_i.size(); ++r){
                charge out_r_charge = right_i[r].first;
                charge out_l_charge = SymmGroup::fuse(out_r_charge, total_delta);             if(!out_left_i.has(out_l_charge)) continue;
                size_t r_size = right_i[r].second;

                size_t o = ret.find_block(out_l_charge, out_r_charge);
                for (size_t w_block = 0; w_block < W.n_blocks(); ++w_block){
                    charge phys_c1 = W.basis().left_charge(w_block);
                    charge phys_c2 = W.basis().right_charge(w_block);

                    charge in_r_charge = SymmGroup::fuse(out_r_charge, -phys_c1);
                    charge in_l_charge = SymmGroup::fuse(in_r_charge, -T_delta);
                    size_t t_block = T.basis().position(in_l_charge, in_r_charge);            if(t_block == T.basis().size()) continue;

                    size_t in_right_offset = in_right_pb(phys_c1, out_r_charge);
                    size_t out_left_offset = out_left_pb(phys_c2, in_l_charge);
                    size_t phys_s1 = W.basis().left_size(w_block);
                    size_t phys_s2 = W.basis().right_size(w_block);
                    Matrix const & wblock = W[w_block];
                    Matrix const & iblock = T[t_block];
                    Matrix & oblock = ret[o];

                    parallel::guard proc(scheduler(o));
                    maquis::dmrg::detail::lb_tensor_mpo(oblock, iblock, wblock,
                                                        out_left_offset, in_right_offset,
                                                        phys_s1, phys_s2, T.basis().left_size(t_block), r_size, access.scale);
                }
            } // right index block
        } // b1
    }

    template<class Matrix, class SymmGroup>
    void rbtm_kernel_allocate(size_t b1,
                              block_matrix<Matrix, SymmGroup> & ret,
                              std::vector<block_matrix<Matrix, SymmGroup> > const & right_mult_mps,
                              MPOTensor<Matrix, SymmGroup> const & mpo,
                              Index<SymmGroup> const & left_i,
                              Index<SymmGroup> const & out_right_i)
    {
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::row_proxy row_proxy;

        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        row_proxy row_b1 = mpo.row(b1);
        for (typename row_proxy::const_iterator row_it = row_b1.begin(); row_it != row_b1.end(); ++row_it) {
            index_type b2 = row_it.index();
            block_matrix<Matrix, SymmGroup> const & T = right_mult_mps[b2];                   if(T.n_blocks() == 0) continue;
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op;                            if(W.n_blocks() == 0) continue;

            charge operator_delta = SymmGroup::fuse(W.basis().right_charge(0), -W.basis().left_charge(0));
            charge        T_delta = SymmGroup::fuse(T.basis().right_charge(0), -T.basis().left_charge(0));
            charge    total_delta = SymmGroup::fuse(operator_delta, -T_delta);

            for(size_t l = 0; l < left_i.size(); ++l){
                charge out_l_charge = left_i[l].first;
                charge out_r_charge = SymmGroup::fuse(out_l_charge, -total_delta);            if(!out_right_i.has(out_r_charge)) continue;
                size_t l_size = left_i[l].second;
                if(ret.find_block(out_l_charge, out_r_charge) == ret.n_blocks())
                    //ret.resize_block(ret.insert_block(Matrix(1,1), out_l_charge, out_r_charge), 
                    //                 l_size, out_right_i.size_of_block(out_r_charge));
                    ret.insert_block(Matrix(l_size, out_right_i.size_of_block(out_r_charge)), out_l_charge, out_r_charge);
            }
        }
        ret.index_sizes();
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void rbtm_kernel_execute(size_t b1,
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
        parallel::scheduler_size_indexed scheduler(ret);

        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::row_proxy row_proxy;

        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        row_proxy row_b1 = mpo.row(b1);
        for (typename row_proxy::const_iterator row_it = row_b1.begin(); row_it != row_b1.end(); ++row_it) {
            index_type b2 = row_it.index();
            block_matrix<Matrix, SymmGroup> const & T = right_mult_mps[b2];                   if(T.n_blocks() == 0) continue;
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op;                            if(W.n_blocks() == 0) continue;

            // charge deltas are constant for all blocks
            charge operator_delta = SymmGroup::fuse(W.basis().right_charge(0), -W.basis().left_charge(0));
            charge        T_delta = SymmGroup::fuse(T.basis().right_charge(0), -T.basis().left_charge(0));
            charge    total_delta = SymmGroup::fuse(operator_delta, -T_delta);

            for (size_t l = 0; l < left_i.size(); ++l){
                charge out_l_charge = left_i[l].first;
                charge out_r_charge = SymmGroup::fuse(out_l_charge, -total_delta);            if(!out_right_i.has(out_r_charge)) continue;
                size_t l_size = left_i[l].second;
                size_t o = ret.find_block(out_l_charge, out_r_charge);

                for (size_t w_block = 0; w_block < W.n_blocks(); ++w_block){
                    charge phys_c1 = W.basis().left_charge(w_block);
                    charge phys_c2 = W.basis().right_charge(w_block);

                    charge in_l_charge = SymmGroup::fuse(out_l_charge, phys_c1); 
                    charge in_r_charge = SymmGroup::fuse(in_l_charge, T_delta); 
                    size_t t_block = T.basis().position(in_l_charge, in_r_charge);            if (t_block == T.basis().size()) continue;

                    size_t in_left_offset = in_left_pb(phys_c1, out_l_charge);
                    size_t out_right_offset = out_right_pb(phys_c2, in_r_charge);
                    size_t phys_s1 = W.basis().left_size(w_block);
                    size_t phys_s2 = W.basis().right_size(w_block);
                    const Matrix & wblock = W[w_block];
                    const Matrix & iblock = T[t_block];
                    Matrix & oblock = ret[o];

                    parallel::guard proc(scheduler(o));
                    maquis::dmrg::detail::rb_tensor_mpo(oblock, iblock, wblock,
                                                        out_right_offset, in_left_offset,
                                                        phys_s1, phys_s2,
                                                        l_size, T.basis().right_size(t_block), access.scale);
                }
            }
        }
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel(size_t b2,
                     ContractionGrid<Matrix, SymmGroup>& contr_grid,
                     Boundary<OtherMatrix, SymmGroup> const & left,
                     std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     DualIndex<SymmGroup> const & ket_basis, // dummy argument for compatiblity with SU2
                     Index<SymmGroup> const & right_i,
                     Index<SymmGroup> const & out_left_i,
                     ProductBasis<SymmGroup> const & in_right_pb,
                     ProductBasis<SymmGroup> const & out_left_pb)
    {
        lbtm_kernel_allocate(b2, contr_grid, left_mult_mps, mpo, right_i, out_left_i);
        lbtm_kernel_execute(b2, contr_grid, left, left_mult_mps, mpo, ket_basis, right_i, out_left_i, in_right_pb, out_left_pb);
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void rbtm_kernel(size_t b1,
                     block_matrix<Matrix, SymmGroup> & ret,
                     Boundary<OtherMatrix, SymmGroup> const & right,
                     std::vector<block_matrix<Matrix, SymmGroup> > const & right_mult_mps,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     DualIndex<SymmGroup> const & ket_basis, // dummy argument for compatiblity with SU2
                     Index<SymmGroup> const & left_i,
                     Index<SymmGroup> const & out_right_i,
                     ProductBasis<SymmGroup> const & in_left_pb,
                     ProductBasis<SymmGroup> const & out_right_pb)
    {
        rbtm_kernel_allocate(b1, ret, right_mult_mps, mpo, left_i, out_right_i);
        rbtm_kernel_execute(b1, ret, right, right_mult_mps, mpo, ket_basis, left_i, out_right_i, in_left_pb, out_right_pb);
    }
    
    } // namespace ablian
} // namespace contraction
#endif
