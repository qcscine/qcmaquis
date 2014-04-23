/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/indexing.h"

#ifdef USE_AMBIENT
#include "dmrg/mp_tensors/impl/ambient.hpp"
#else
#include "dmrg/mp_tensors/impl/alps.hpp"
#endif

namespace contraction {

    // output/input: left_i for bra_tensor, right_i for ket_tensor
    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<OtherMatrix, SymmGroup>
    overlap_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                      MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                      block_matrix<OtherMatrix, SymmGroup> const & left,
                      block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
    {
        if (localop != NULL)
            throw std::runtime_error("Not implemented!");
        
        assert(ket_tensor.phys_i == bra_tensor.phys_i);
        
        bra_tensor.make_left_paired();
        
        block_matrix<OtherMatrix, SymmGroup> t1;
        block_matrix<Matrix, SymmGroup> t3;
        ket_tensor.make_right_paired();
        gemm(left, ket_tensor.data(), t1);
        
        reshape_right_to_left_new(ket_tensor.site_dim(), bra_tensor.row_dim(), ket_tensor.col_dim(),
                                  t1, t3);
        gemm(transpose(conjugate(bra_tensor.data())), t3, t1);
        return t1;

        // original:
        // t3 = transpose(t3);
        // gemm(t3, t2, t1);
        // return transpose(t1);
    }
    
    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<OtherMatrix, SymmGroup>
    overlap_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                       MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                       block_matrix<OtherMatrix, SymmGroup> const & right,
                       block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
    {
        if (localop != NULL)
            throw std::runtime_error("Not implemented!");
        
        assert(ket_tensor.phys_i == bra_tensor.phys_i);
        
        bra_tensor.make_right_paired();
        ket_tensor.make_left_paired();
        
        block_matrix<OtherMatrix, SymmGroup> t1;
        block_matrix<Matrix, SymmGroup> t3;
        gemm(ket_tensor.data(), transpose(right), t1);
        reshape_left_to_right_new(ket_tensor.site_dim(), ket_tensor.row_dim(), bra_tensor.col_dim(), t1, t3);
        gemm(conjugate(bra_tensor.data()), transpose(t3), t1);

        return t1;
    }
    
    // SK: New version which generates same output but uses right-paired input.
    //     The charge delta optimization is indepent from the changes needed to
    //     skip the preceding reshapes.
    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel(size_t b2,
                     ContractionGrid<Matrix, SymmGroup>& contr_grid,
                     Boundary<OtherMatrix, SymmGroup> const & left,
                     std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     Index<SymmGroup> const & physical_i,
                     Index<SymmGroup> const & right_i,
                     Index<SymmGroup> const & out_left_i,
                     ProductBasis<SymmGroup> const & in_right_pb,
                     ProductBasis<SymmGroup> const & out_left_pb)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;

        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();

            block_matrix<Matrix, SymmGroup> const & T = left_mult_mps[b1];
            if (T.n_blocks() == 0) continue;
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            block_matrix<Matrix, SymmGroup> const & W = access.op;
            if (W.n_blocks() == 0) continue;

            // charge deltas are constant for all blocks
            charge operator_delta = SymmGroup::fuse(W.right_basis()[0].first, -W.left_basis()[0].first);
            charge        T_delta = SymmGroup::fuse(T.right_basis()[0].first, -T.left_basis()[0].first);
            charge    total_delta = SymmGroup::fuse(operator_delta, -T_delta);
        
            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);

            for (size_t r = 0; r < right_i.size(); ++r)
            {
                charge out_r_charge = right_i[r].first;
                charge out_l_charge = SymmGroup::fuse(out_r_charge, total_delta);
                size_t r_size = right_i[r].second;

                if (!out_left_i.has(out_l_charge)) continue;

                size_t o = ret.find_block(out_l_charge, out_r_charge);
                if ( o == ret.n_blocks() ) {
                    o = ret.insert_block(Matrix(1,1), out_l_charge, out_r_charge);
                    ret.resize_block(out_l_charge, out_r_charge, out_left_i.size_of_block(out_l_charge), r_size);
                }

                for (size_t w_block = 0; w_block < W.n_blocks(); ++w_block)
                {
                    charge phys_c1 = W.left_basis()[w_block].first;
                    charge phys_c2 = W.right_basis()[w_block].first;

                    charge in_r_charge = SymmGroup::fuse(out_r_charge, -phys_c1);
                    size_t t_block = T.right_basis().position(in_r_charge);
                    if (t_block == T.right_basis().size()) continue;
 
                    charge in_l_charge = T.left_basis()[t_block].first;

                    size_t in_right_offset = in_right_pb(phys_c1, out_r_charge);
                    size_t out_left_offset = out_left_pb(phys_c2, in_l_charge);

                    size_t phys_s1 = W.left_basis()[w_block].second;
                    size_t phys_s2 = W.right_basis()[w_block].second;
                    Matrix const & wblock = W[w_block];
                    Matrix const & iblock = T[t_block];
                    Matrix & oblock = ret[o];

                    maquis::dmrg::detail::lb_tensor_mpo(oblock, iblock, wblock,
                            out_left_offset, in_right_offset,
                            phys_s1, phys_s2, T.left_basis()[t_block].second, r_size, access.scale);
                }
            } // right index block
        } // b1
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    rbtm_kernel(size_t b1,
                Boundary<OtherMatrix, SymmGroup> const & right,
                std::vector<block_matrix<Matrix, SymmGroup> > const & right_mult_mps,
                MPOTensor<Matrix, SymmGroup> const & mpo,
                Index<SymmGroup> const & physical_i,
                Index<SymmGroup> const & left_i,
                Index<SymmGroup> const & right_i,
                Index<SymmGroup> const & out_right_i,
                ProductBasis<SymmGroup> const & in_left_pb,
                ProductBasis<SymmGroup> const & out_right_pb)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;

        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        block_matrix<Matrix, SymmGroup> ret;

        row_proxy row_b1 = mpo.row(b1);
        for (typename row_proxy::const_iterator row_it = row_b1.begin(); row_it != row_b1.end(); ++row_it) {
            index_type b2 = row_it.index();

            block_matrix<Matrix, SymmGroup> const & T = right_mult_mps[b2];
            if (T.n_blocks() == 0) continue;
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            block_matrix<Matrix, SymmGroup> const & W = access.op;
            if (W.n_blocks() == 0) continue;

            // charge deltas are constant for all blocks
            charge operator_delta = SymmGroup::fuse(W.right_basis()[0].first, -W.left_basis()[0].first);
            charge        T_delta = SymmGroup::fuse(T.right_basis()[0].first, -T.left_basis()[0].first);
            charge    total_delta = SymmGroup::fuse(operator_delta, -T_delta);

            for (size_t l = 0; l < left_i.size(); ++l)
            {
                charge out_l_charge = left_i[l].first;
                size_t l_size = left_i[l].second;
                charge out_r_charge = SymmGroup::fuse(out_l_charge, -total_delta);

                if (!out_right_i.has(out_r_charge)) continue;

                size_t o = ret.find_block(out_l_charge, out_r_charge);
                if ( o == ret.n_blocks() ) {
                    o = ret.insert_block(Matrix(1,1), out_l_charge, out_r_charge);
                    ret.resize_block(out_l_charge, out_r_charge, l_size, out_right_i.size_of_block(out_r_charge));
                }

                for (size_t w_block = 0; w_block < W.n_blocks(); ++w_block)
                {
                    charge phys_c1 = W.left_basis()[w_block].first;
                    charge phys_c2 = W.right_basis()[w_block].first;

                    charge in_l_charge = SymmGroup::fuse(out_l_charge,  phys_c1); 
                    size_t t_block = T.left_basis().position(in_l_charge);
                    if (t_block == T.left_basis().size()) continue;

                    charge in_r_charge = T.right_basis()[t_block].first;
                    assert(right_i.has(in_r_charge));

                    size_t in_left_offset = in_left_pb(phys_c1, out_l_charge);
                    size_t out_right_offset = out_right_pb(phys_c2, in_r_charge);
                    size_t phys_s1 = W.left_basis()[w_block].second;
                    size_t phys_s2 = W.right_basis()[w_block].second;
                    const Matrix & wblock = W[w_block];
                    const Matrix & iblock = T[t_block];
                    Matrix & oblock = ret[o];

                    maquis::dmrg::detail::rb_tensor_mpo(oblock, iblock, wblock,
                            out_right_offset, in_left_offset,
                            phys_s1, phys_s2,
                            l_size, T.right_basis()[t_block].second, access.scale);

                }
            }
        }
        return ret;
    }

    // note: this function changes the internal structure of Boundary,
    //       each block is transposed
    template<class Matrix, class OtherMatrix, class SymmGroup>
    Boundary<Matrix, SymmGroup>
    left_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                             Boundary<OtherMatrix, SymmGroup> const & left,
                             MPOTensor<Matrix, SymmGroup> const & mpo,
                             Index<SymmGroup> const * in_low = NULL)
    {
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        if (in_low == NULL)
            in_low = &mps.row_dim();
        
        std::vector<block_matrix<Matrix, SymmGroup> > t(left.aux_dim());
        int loop_max = left.aux_dim();

        {
            select_proc(ambient::scope_t::common);
            mps.make_right_paired();
            storage::hint(mps);
            gemm_trim_left(transpose(left[0]), mps.data(), t[0]);
        }
        parallel_for(int b, range(1,loop_max), {
            select_proc(ambient::scope::permute(b,mpo.placement_l)); 
            gemm_trim_left(transpose(left[b]), mps.data(), t[b]);
        });

        Index<SymmGroup> physical_i = mps.site_dim(), left_i = *in_low, right_i = mps.col_dim(),
                                      out_left_i = physical_i * left_i;
        ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> in_right_pb(physical_i, right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));
        
        loop_max = mpo.col_dim();

        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, left.aux_dim(), mpo.col_dim());

        parallel_for(int b2, range(0,loop_max), {
            select_proc(ambient::scope::permute(b2,mpo.placement_r));
            lbtm_kernel(b2, contr_grid, left, t, mpo, physical_i, right_i, out_left_i, in_right_pb, out_left_pb);
        });

        return contr_grid.make_boundary();
    }
    
    template<class Matrix, class OtherMatrix, class SymmGroup>
    Boundary<Matrix, SymmGroup>
    right_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                              Boundary<OtherMatrix, SymmGroup> const & right,
                              MPOTensor<Matrix, SymmGroup> const & mpo,
                              Index<SymmGroup> const * in_low = NULL)
    {
        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        if (in_low == NULL)
            in_low = &mps.col_dim();
        
        std::vector<block_matrix<Matrix, SymmGroup> > t(right.aux_dim());
        int loop_max = right.aux_dim();
        
        {
            select_proc(storage::scope_t::common);
            mps.make_left_paired();
            storage::hint(mps);
        }
        parallel_for(int b, range(0,loop_max), {
            select_proc(ambient::scope::permute(b,mpo.placement_r));
            gemm_trim_right(mps.data(), right[b], t[b]);
        });
        
        Index<SymmGroup> physical_i = mps.site_dim(), left_i = mps.row_dim(), right_i = *in_low,
                         out_right_i = adjoin(physical_i) * right_i;

        ProductBasis<SymmGroup> in_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> out_right_pb(physical_i, right_i,
                                             boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                 -boost::lambda::_1, boost::lambda::_2));
        Boundary<Matrix, SymmGroup> ret;
        ret.resize(mpo.row_dim());
        
        loop_max = mpo.row_dim();

        omp_for(int b1 = 0; b1 < loop_max; ++b1) {
            select_proc(ambient::scope::permute(b1,mpo.placement_l));
            ret[b1] = rbtm_kernel(b1, right, t, mpo, physical_i, left_i, right_i, out_right_i, in_left_pb, out_right_pb);
        }

        return ret;
    }
    
    template<class Matrix, class OtherMatrix, class SymmGroup>
    Boundary<OtherMatrix, SymmGroup>
    overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          Boundary<OtherMatrix, SymmGroup> const & left,
                          MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        #ifdef AMBIENT_TRACKING
        ambient::overseer::log::region("parallel::overlap_mpo_left_step");
        #endif

        typedef typename SymmGroup::charge charge;

        std::vector<block_matrix<Matrix, SymmGroup> > t(left.aux_dim());
        {
            // Make a copy of ket_tensor to avoid reshaping back to left
            MPSTensor<Matrix, SymmGroup> ket_cpy = ket_tensor;
            int loop_max = left.aux_dim();
            {
                select_proc(ambient::scope_t::common);
                ket_cpy.make_right_paired();
                storage::hint(ket_cpy);
                gemm_trim_left(transpose(left[0]), ket_cpy.data(), t[0]);
            }
            parallel_for(int b, range(1,loop_max), {
                select_proc(ambient::scope::permute(b,mpo.placement_l));
                gemm_trim_left(transpose(left[b]), ket_cpy.data(), t[b]);
            });
        }

        Index<SymmGroup> const & left_i = bra_tensor.row_dim();
        Index<SymmGroup> const & right_i = ket_tensor.col_dim();
        Index<SymmGroup> out_left_i = ket_tensor.site_dim() * left_i;
        ProductBasis<SymmGroup> out_left_pb(ket_tensor.site_dim(), left_i);
        ProductBasis<SymmGroup> in_right_pb(ket_tensor.site_dim(), right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));

        int loop_max = mpo.col_dim();

        bra_tensor.make_left_paired();
        block_matrix<Matrix, SymmGroup> bra_conj = conjugate(bra_tensor.data());

        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, left.aux_dim(), mpo.col_dim());
        parallel_for(int b2, range(0,loop_max), {
            select_proc(ambient::scope::permute(b2,mpo.placement_r));
            lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.site_dim(), right_i, out_left_i, in_right_pb, out_left_pb);
            contr_grid.multiply_column_trans(b2, bra_conj);
        });
        #ifdef AMBIENT_TRACKING
        ambient::overseer::log::region("serial::continue");
        #endif

        return contr_grid.make_boundary();
    }
    
    template<class Matrix, class OtherMatrix, class SymmGroup>
    Boundary<OtherMatrix, SymmGroup>
    overlap_mpo_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                           MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                           Boundary<OtherMatrix, SymmGroup> const & right,
                           MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        typedef typename SymmGroup::charge charge;

        #ifdef AMBIENT_TRACKING
        ambient::overseer::log::region("parallel::overlap_mpo_right_step");
        #endif

        std::vector<block_matrix<Matrix, SymmGroup> > t(right.aux_dim());
        {
            // Make a copy of ket_tensor to avoid reshaping back to right
            MPSTensor<Matrix, SymmGroup> ket_cpy = ket_tensor;
            ket_cpy.make_left_paired();
            int loop_max = right.aux_dim();

            storage::hint(ket_cpy, storage::scope_t::common);
            parallel_for(int b, range(0,loop_max), {
                select_proc(ambient::scope::permute(b,mpo.placement_r));
                gemm_trim_right(ket_cpy.data(), right[b], t[b]);
            });
        }

        Index<SymmGroup> const & left_i = ket_tensor.row_dim();
        Index<SymmGroup> const & right_i = bra_tensor.col_dim();
        Index<SymmGroup> out_right_i = adjoin(ket_tensor.site_dim()) * right_i;
        ProductBasis<SymmGroup> in_left_pb(ket_tensor.site_dim(), left_i);
        ProductBasis<SymmGroup> out_right_pb(ket_tensor.site_dim(), right_i,
                                             boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                 -boost::lambda::_1, boost::lambda::_2));
        Boundary<Matrix, SymmGroup> ret;
        ret.resize(mpo.row_dim());

        //ket_tensor.make_right_paired();
        std::size_t loop_max = mpo.row_dim();

        bra_tensor.make_right_paired();
        block_matrix<Matrix, SymmGroup> bra_conj = conjugate(bra_tensor.data());
        omp_for(size_t b1 = 0; b1 < loop_max; ++b1) {
            select_proc(ambient::scope::permute(b1,mpo.placement_l));
            block_matrix<Matrix, SymmGroup> tmp;
            tmp = rbtm_kernel(b1, right, t, mpo, ket_tensor.site_dim(), left_i, right_i, out_right_i, in_left_pb, out_right_pb);
            gemm(tmp, transpose(bra_conj), ret[b1]);
        }
        #ifdef AMBIENT_TRACKING
        ambient::overseer::log::region("serial::continue");
        #endif

        return ret;
    }
    
    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
                Boundary<OtherMatrix, SymmGroup> const & left,
                Boundary<OtherMatrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        typedef typename SymmGroup::charge charge;

        ket_tensor.make_right_paired();
        
        std::vector<block_matrix<Matrix, SymmGroup> > t(left.aux_dim());
        int loop_max = left.aux_dim();
        {
            select_proc(storage::scope_t::common);
            storage::hint(ket_tensor);
            storage::migrate(left[0]);
            gemm_trim_left(transpose(left[0]), ket_tensor.data(), t[0]);
        }
        parallel_for(int b1, range(1,loop_max), {
            select_proc(ambient::scope::permute(b1,mpo.placement_l));
            gemm_trim_left(transpose(left[b1]), ket_tensor.data(), t[b1]);
        });

        Index<SymmGroup> const & physical_i = ket_tensor.site_dim(),
                               & left_i = ket_tensor.row_dim(),
                               & right_i = ket_tensor.col_dim(),
                                 out_left_i = physical_i * left_i;
        ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> in_right_pb(physical_i, right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));

        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, left.aux_dim(), mpo.col_dim());

        loop_max = mpo.col_dim();
        parallel_for(int b2, range(0,loop_max), {
            select_proc(ambient::scope::permute(b2,mpo.placement_r));
            lbtm_kernel(b2, contr_grid, left, t, mpo, physical_i, right_i, out_left_i, in_right_pb, out_left_pb);
            contr_grid.multiply_column(b2, right[b2]);
        });
         
        MPSTensor<Matrix, SymmGroup> ret;
        ret.phys_i = ket_tensor.site_dim();
        ret.left_i = ket_tensor.row_dim();
        ret.right_i = ket_tensor.col_dim();
        ret.data() = contr_grid.reduce();
        return ret;
    }
    
    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    site_ortho_boundaries(MPSTensor<Matrix, SymmGroup> const & mps,
                          MPSTensor<Matrix, SymmGroup> const & ortho_mps,
                          block_matrix<OtherMatrix, SymmGroup> const & ortho_left,
                          block_matrix<OtherMatrix, SymmGroup> const & ortho_right)
    {
        ortho_mps.make_right_paired();
        block_matrix<Matrix, SymmGroup> t, t2, t3;
        gemm(ortho_left, ortho_mps.data(), t);
        reshape_right_to_left_new(mps.site_dim(),
                                  ortho_left.left_basis(), ortho_mps.col_dim(),
                                  t, t2);
        gemm(t2, transpose(ortho_right), t3);
        
        mps.make_left_paired();
        t = mps.data();
        reshape_and_pad_left(mps.site_dim(),
                             ortho_left.left_basis(), ortho_right.left_basis(),
                             mps.row_dim(), mps.col_dim(),
                             t3, t);
        
        MPSTensor<Matrix, SymmGroup> t4(mps.site_dim(),
                                        mps.row_dim(), mps.col_dim(),
                                        t, LeftPaired);
        return t4;
    }
    
    template<class Matrix, class OtherMatrix, class SymmGroup>
    std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
    predict_new_state_l2r_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                MPOTensor<Matrix, SymmGroup> const & mpo,
                                Boundary<OtherMatrix, SymmGroup> const & left,
                                Boundary<OtherMatrix, SymmGroup> const & right,
                                double alpha, double cutoff, std::size_t Mmax)
    {
        mps.make_left_paired();
        block_matrix<Matrix, SymmGroup> dm;
        gemm(mps.data(), transpose(conjugate(mps.data())), dm);
        
        Boundary<Matrix, SymmGroup> half_dm = left_boundary_tensor_mpo(mps, left, mpo);
        
        mps.make_left_paired();
        for (std::size_t b = 0; b < half_dm.aux_dim(); ++b)
        {
            block_matrix<Matrix, SymmGroup> tdm;
            gemm(half_dm[b], transpose(conjugate(half_dm[b])), tdm);
            
            
            tdm *= alpha;
            for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                if (mps.data().left_basis().has(tdm.left_basis()[k].first))
                    dm.match_and_add_block(tdm[k],
                                           tdm.left_basis()[k].first,
                                           tdm.right_basis()[k].first);
            }
        }
        mps.make_left_paired();
        assert( weak_equal(dm.left_basis(), mps.data().left_basis()) );
        
        block_matrix<Matrix, SymmGroup> U;
        block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
        truncation_results trunc = heev_truncate(dm, U, S, cutoff, Mmax);
      
        MPSTensor<Matrix, SymmGroup> ret = mps;
        ret.replace_left_paired(U);
        return std::make_pair(ret, trunc);
    }
    
    template<class Matrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    predict_lanczos_l2r_sweep(MPSTensor<Matrix, SymmGroup> B,
                              MPSTensor<Matrix, SymmGroup> const & psi,
                              MPSTensor<Matrix, SymmGroup> const & A)
    {
        psi.make_left_paired();
        A.make_left_paired();
        
        block_matrix<Matrix, SymmGroup> tmp;
        gemm(transpose(conjugate(A.data())), psi.data(), tmp);
        B.multiply_from_left(tmp);
        
        return B;
    }
    
    template<class Matrix, class OtherMatrix, class SymmGroup>
    std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
    predict_new_state_r2l_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                MPOTensor<Matrix, SymmGroup> const & mpo,
                                Boundary<OtherMatrix, SymmGroup> const & left,
                                Boundary<OtherMatrix, SymmGroup> const & right,
                                double alpha, double cutoff, std::size_t Mmax)
    {
        mps.make_right_paired();
        block_matrix<Matrix, SymmGroup> dm;
        gemm(transpose(conjugate(mps.data())), mps.data(), dm);
            
        Boundary<Matrix, SymmGroup> half_dm = right_boundary_tensor_mpo(mps, right, mpo);
        
        mps.make_right_paired();
        for (std::size_t b = 0; b < half_dm.aux_dim(); ++b)
        {
            block_matrix<Matrix, SymmGroup> tdm;
            gemm(transpose(conjugate(half_dm[b])), half_dm[b], tdm);
            
            tdm *= alpha;
            for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                if (mps.data().right_basis().has(tdm.right_basis()[k].first))
                    dm.match_and_add_block(tdm[k],
                                           tdm.left_basis()[k].first,
                                           tdm.right_basis()[k].first);
            }
        }
        
        mps.make_right_paired();
        assert( weak_equal(dm.right_basis(), mps.data().right_basis()) );
        
        block_matrix<Matrix, SymmGroup> U;
        block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
        truncation_results trunc = heev_truncate(dm, U, S, cutoff, Mmax);
        
        MPSTensor<Matrix, SymmGroup> ret = mps;
        ret.replace_right_paired(adjoint(U));
        return std::make_pair(ret, trunc);
    }
    
    template<class Matrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    predict_lanczos_r2l_sweep(MPSTensor<Matrix, SymmGroup> B,
                              MPSTensor<Matrix, SymmGroup> const & psi,
                              MPSTensor<Matrix, SymmGroup> const & A)
    {
        psi.make_right_paired();
        A.make_right_paired();
        
        block_matrix<Matrix, SymmGroup> tmp;
        gemm(psi.data(), transpose(conjugate(A.data())), tmp);
        
        B.multiply_from_right(tmp);
        
        return B;
    }
    
    template<class Matrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    multiply_with_twosite(block_matrix<Matrix, SymmGroup> const & vec,
                          block_matrix<Matrix, SymmGroup> const & op,
                          Index<SymmGroup> const & left_i,
                          Index<SymmGroup> const & right_i,
                          Index<SymmGroup> const & phys_i)
    {
        typedef typename SymmGroup::charge charge;
        
        Index<SymmGroup> l_index = phys_i * left_i, r_index = adjoin(phys_i) * right_i;
        common_subset(l_index, r_index);
        block_matrix<Matrix, SymmGroup> ret;
        
        ProductBasis<SymmGroup> left_pb(phys_i, left_i);
        ProductBasis<SymmGroup> right_pb(phys_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
        ProductBasis<SymmGroup> phys_pb(phys_i, phys_i);
        
        for (size_t ls = 0; ls < left_i.size(); ++ls)
            for (size_t rs = 0; rs < right_i.size(); ++rs)
                for (size_t lps = 0; lps < phys_i.size(); ++lps)
                    for (size_t rps = 0; rps < phys_i.size(); ++rps)
                        for (size_t ilps = 0; ilps < phys_i.size(); ++ilps)
                            for (size_t irps = 0; irps < phys_i.size(); ++irps)
                            {
                                charge lc = left_i[ls].first, rc = right_i[rs].first;
                                charge lpc = phys_i[lps].first, rpc = phys_i[rps].first;
                                charge ilpc = phys_i[ilps].first, irpc = phys_i[irps].first;
                                
                                charge left_vec_charge = SymmGroup::fuse(ilpc, lc);
                                charge right_vec_charge = SymmGroup::fuse(-irpc, rc);
                                
                                charge left_out_charge = SymmGroup::fuse(lpc, lc);
                                charge right_out_charge = SymmGroup::fuse(-rpc, rc);
                                
                                if (left_out_charge != right_out_charge)
                                    continue;
                                if (left_vec_charge != right_vec_charge)
                                    continue;
                                if (! vec.has_block(left_vec_charge, right_vec_charge) )
                                    continue;
                                if (SymmGroup::fuse(lpc, rpc) != SymmGroup::fuse(ilpc, irpc))
                                    continue;

                                if (! ret.has_block(left_out_charge, right_out_charge))
                                    ret.insert_block(Matrix(l_index.size_of_block(left_out_charge), r_index.size_of_block(right_out_charge), 0),
                                                     left_out_charge, right_out_charge);

                                charge both_charge = SymmGroup::fuse(lpc, rpc);
                                
                                assert( op.has_block(both_charge, both_charge) );
                                
                                maquis::dmrg::detail::mwt(ret(left_out_charge, right_out_charge), 
                                                          vec(left_vec_charge, right_vec_charge), 
                                                          op(both_charge, both_charge), 
                                                          left_pb(lpc, lc), right_pb(rpc, rc),
                                                          left_pb(ilpc, lc), right_pb(irpc, rc),
                                                          phys_pb(ilpc, irpc), phys_pb(lpc, rpc),
                                                          left_i[ls].second, right_i[rs].second,
                                                          phys_i[lps].second, phys_i[rps].second,
                                                          phys_i[ilps].second, phys_i[irps].second);
                            }
        
        return ret;
    }
    
    template<class Matrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    multiply_with_op(MPSTensor<Matrix, SymmGroup> const & mps, block_matrix<Matrix, SymmGroup> const & op)
    {
        typedef typename SymmGroup::charge charge;
        
        mps.make_left_paired();
        block_matrix<Matrix, SymmGroup> const & vec = mps.data();
        
        Index<SymmGroup> const& phys_i  = mps.site_dim();
        Index<SymmGroup> const& left_i  = mps.row_dim();
        Index<SymmGroup> const& right_i = mps.col_dim();
        
        ProductBasis<SymmGroup> left_pb(phys_i, left_i);
        
        block_matrix<Matrix, SymmGroup> ret_vec;
        Index<SymmGroup> out_left_i;

        for (size_t ls = 0; ls < left_i.size(); ++ls)
            for (size_t rs = 0; rs < right_i.size(); ++rs)
                for (size_t s1 = 0; s1 < phys_i.size(); ++s1)
                    for (size_t s2 = 0; s2 < phys_i.size(); ++s2) {
                        
                        charge lc  = left_i[ls].first, rc  = right_i[rs].first;
                        charge s1c = phys_i[s1].first, s2c = phys_i[s2].first;
                        
                        if (! op.has_block(s1c, s2c) )
                            continue;

                        charge phys_diff = SymmGroup::fuse(s2c, -s1c);
                        
                        charge left_vec_charge  = SymmGroup::fuse(s1c, lc);
                        charge right_vec_charge = rc;

                        if (! vec.has_block(left_vec_charge, right_vec_charge) )
                            continue;
                        
                        if (! out_left_i.has(lc))
                            out_left_i.insert(left_i[ls]);
                        
                        charge left_out_charge  = SymmGroup::fuse(s2c, lc);
                        charge right_out_charge = SymmGroup::fuse(rc, phys_diff);
                        
                        if (! ret_vec.has_block(left_out_charge, right_out_charge) )
                            ret_vec.insert_block(new Matrix(left_pb.size(s2c, lc), right_i[rs].second, 0), left_out_charge, right_out_charge);
                        
                        Matrix & oblock         = ret_vec(left_out_charge, right_out_charge);
                        Matrix const & iblock   = vec(left_vec_charge, right_vec_charge);
                        Matrix const & op_block = op(s1c, s2c);
                        
                        size_t i_l_offset = left_pb(s1c, lc);
                        size_t i_r_offset = 0;
                        size_t l_offset = left_pb(s2c, lc);
                        size_t r_offset = 0;
                        size_t i_op_offset = 0;
                        size_t op_offset = 0;
                        
#ifdef USE_AMBIENT
                        printf("UNOPTIMIZED FUNCTION (MULTIPLY WITH OP!)\n");
#endif
                        for (size_t ss1 = 0; ss1 < phys_i[s1].second; ++ss1) {
                            for (size_t ss2 = 0; ss2 < phys_i[s2].second; ++ss2) {
                                size_t o_l_start = l_offset + ss2*left_i[ls].second;
                                size_t o_r_start = r_offset;
                                
                                size_t i_l_start = i_l_offset + ss1*left_i[ls].second;
                                size_t i_r_start = i_r_offset;
                                
                                typename Matrix::value_type const & op_val = op_block(i_op_offset + ss1,  op_offset + ss2);
                                
                                // TODO: replace by kernel
                                for (size_t rr = 0; rr < right_i[rs].second; ++rr) {
                                    for (size_t ll = 0; ll < left_i[ls].second; ++ll) {
                                        oblock(o_l_start+ll, o_r_start+rr) += iblock(i_l_start+ll, i_r_start+rr) * op_val;
                                    }
                                }
                            }
                        }
                    }
        
        
        MPSTensor<Matrix, SymmGroup> ret(phys_i, out_left_i, ret_vec.right_basis(), ret_vec, LeftPaired);
        assert( ret.reasonable() );
        return ret;
    }
    
    // tested only for diagonal operator
    template<class Matrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    local_op(MPSTensor<Matrix, SymmGroup> const & mps,
             block_matrix<Matrix, SymmGroup> const & op)
    {
        typedef typename SymmGroup::charge charge;
        
        mps.make_left_paired();
        block_matrix<Matrix, SymmGroup> const & vec = mps.data();
        
        Index<SymmGroup> phys_i = mps.site_dim();
        Index<SymmGroup> left_i = mps.row_dim();
        Index<SymmGroup> right_i = mps.col_dim();
        
        ProductBasis<SymmGroup> left_pb(phys_i, left_i);
        
        MPSTensor<Matrix, SymmGroup> ret(phys_i, left_i, right_i, false);
        block_matrix<Matrix, SymmGroup> & ret_vec = ret.data();
        
        for (size_t ls = 0; ls < left_i.size(); ++ls)
            for (size_t rs = 0; rs < right_i.size(); ++rs)
                for (size_t s1 = 0; s1 < phys_i.size(); ++s1)
                    for (size_t s2 = 0; s2 < phys_i.size(); ++s2) {
                        
                        charge lc = left_i[ls].first, rc = right_i[rs].first;
                        charge s1c = phys_i[s1].first, s2c = phys_i[s2].first;
                        
                        charge left_vec_charge = SymmGroup::fuse(s1c, lc);
                        charge right_vec_charge = rc;
                        charge left_out_charge = SymmGroup::fuse(s2c, lc);
                        charge right_out_charge = right_vec_charge;
                        
                        if (! vec.has_block(left_vec_charge, right_vec_charge) )
                            continue;
                        if (! op.has_block(s1c, s2c) )
                            continue;
                        if (! ret_vec.has_block(left_out_charge, right_out_charge) )
                            ret_vec.insert_block(new Matrix(left_pb.size(s2c, lc), right_i[rs].second, 0), left_out_charge, right_out_charge);
                        
                        Matrix & oblock = ret_vec(left_out_charge, right_out_charge);
                        Matrix const & iblock = vec(left_vec_charge, right_vec_charge);
                        Matrix const & op_block = op(s1c, s2c);
                        
                        size_t i_l_offset = left_pb(s1c, lc);
                        size_t i_r_offset = 0;
                        size_t l_offset = left_pb(s2c, lc);
                        size_t r_offset = 0;
                        size_t i_op_offset = 0;
                        size_t op_offset = 0;
                        
                        for (size_t ll = 0; ll < left_i[ls].second; ++ll)
                            for (size_t rr = 0; rr < right_i[rs].second; ++rr)
                                for (size_t ss1 = 0; ss1 < phys_i[s1].second; ++ss1)
                                    for (size_t ss2 = 0; ss2 < phys_i[s2].second; ++ss2) {
#ifndef NDEBUG
                                        oblock(l_offset + ss2*left_i[ls].second + ll,
                                               r_offset + rr);
                                        iblock(i_l_offset + ss1*left_i[ls].second + ll,
                                               i_r_offset + rr);
                                        op_block(i_op_offset + ss1,
                                                 op_offset + ss2);
#endif
                                        oblock(l_offset + ss2*left_i[ls].second + ll,
                                               r_offset + rr) += 
                                        iblock(i_l_offset + ss1*left_i[ls].second + ll,
                                               i_r_offset + rr) *
                                        op_block(i_op_offset + ss1,
                                                 op_offset + ss2);
                                        
                                    }
                    }
        
        
        assert( ret.reasonable() );
        return ret;
    }

    template<class Matrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    density_matrix(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                   MPSTensor<Matrix, SymmGroup> const & ket_tensor)
    {
        typedef typename SymmGroup::charge charge;
        
        assert( bra_tensor.row_dim() == ket_tensor.row_dim() );
        assert( bra_tensor.col_dim() == ket_tensor.col_dim() );
        assert( bra_tensor.site_dim() == ket_tensor.site_dim() );
        
        bra_tensor.make_left_paired();
        ket_tensor.make_left_paired();
        
        Index<SymmGroup> phys_i = ket_tensor.site_dim();
        Index<SymmGroup> left_i = ket_tensor.row_dim();
        Index<SymmGroup> right_i = ket_tensor.col_dim();
        
        ProductBasis<SymmGroup> left_pb(phys_i, left_i);
        
        block_matrix<Matrix, SymmGroup> const & bra_vec = bra_tensor.data();
        block_matrix<Matrix, SymmGroup> const & ket_vec = ket_tensor.data();
        
        block_matrix<Matrix, SymmGroup> ret(phys_i, phys_i);
        
        for (size_t ls = 0; ls < left_i.size(); ++ls)
            for (size_t rs = 0; rs < right_i.size(); ++rs)
                for (size_t s1 = 0; s1 < phys_i.size(); ++s1)
                    for (size_t s2 = 0; s2 < phys_i.size(); ++s2) {
                        
                        charge lc = left_i[ls].first, rc = right_i[rs].first;
                        charge s1c = phys_i[s1].first, s2c = phys_i[s2].first;
                        
                        charge left_ket_charge = SymmGroup::fuse(s1c, lc);
                        charge left_bra_charge = SymmGroup::fuse(s2c, lc);
                        charge right_charge = rc;
                        
                        if (! ket_vec.has_block(left_ket_charge, right_charge) )
                            continue;
                        if (! bra_vec.has_block(left_bra_charge, right_charge) )
                            continue;
                        if (! ret.has_block(s1c, s2c) )
                            ret.insert_block(new Matrix(phys_i[s1].second, phys_i[s2].second, 0), s1c, s2c);
                        
                        Matrix & oblock = ret(s1c, s2c);
                        Matrix const & ket_block = ket_vec(left_ket_charge, right_charge);
                        Matrix const & bra_block = bra_vec(left_bra_charge, right_charge);
                        
                        size_t l_ket_offset = left_pb(s1c, lc);
                        size_t l_bra_offset = left_pb(s2c, lc);
                        
                        for (size_t ll = 0; ll < left_i[ls].second; ++ll)
                            for (size_t rr = 0; rr < right_i[rs].second; ++rr)
                                for (size_t ss1 = 0; ss1 < phys_i[s1].second; ++ss1)
                                    for (size_t ss2 = 0; ss2 < phys_i[s2].second; ++ss2) {
                                        
                                        oblock(ss1, ss2) +=
                                        utils::conj(bra_block(l_bra_offset + ss2*left_i[ls].second + ll, rr)) *
                                        ket_block(l_ket_offset + ss1*left_i[ls].second + ll, rr);
                                        
                                    }
                    }
        return ret;
    }

    template<class Matrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    density_matrix_2(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                     MPSTensor<Matrix, SymmGroup> const & ket_tensor)
    {
        typedef typename SymmGroup::charge charge;
        
        assert( bra_tensor.row_dim() == ket_tensor.row_dim() );
        assert( bra_tensor.col_dim() == ket_tensor.col_dim() );
        assert( bra_tensor.site_dim() == ket_tensor.site_dim() );
        
        bra_tensor.make_left_paired();
        ket_tensor.make_left_paired();
        
        Index<SymmGroup> phys_i = ket_tensor.site_dim();
        Index<SymmGroup> left_i = ket_tensor.row_dim();
        Index<SymmGroup> right_i = ket_tensor.col_dim();
        
        block_matrix<Matrix, SymmGroup> ket_mat = reshape_left_to_physleft(phys_i, left_i, right_i, ket_tensor.data());
        block_matrix<Matrix, SymmGroup> bra_mat = reshape_left_to_physleft(phys_i, left_i, right_i, bra_tensor.data());

        block_matrix<Matrix, SymmGroup> dm;
        gemm(ket_mat, transpose(conjugate(bra_mat)), dm);
        
        return dm;
    }

}

#endif
