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

#ifndef CONTRACTIONS_SU2_APPLY_OP_HPP
#define CONTRACTIONS_SU2_APPLY_OP_HPP

#include "dmrg/block_matrix/symmetry/gsl_coupling.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/contractions/non-abelian/functors.h"
#include "dmrg/mp_tensors/contractions/non-abelian/micro_kernels.hpp"

#include "dmrg/mp_tensors/contractions/non-abelian/gemm.hpp"

namespace contraction {
namespace SU2 {

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel(size_t b2,
                     ContractionGrid<Matrix, SymmGroup>& contr_grid,
                     Boundary<OtherMatrix, SymmGroup> const & left,
                     BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> const & left_mult_mps,
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

            block_matrix<Matrix, SymmGroup> local;
            block_matrix<Matrix, SymmGroup> const & T = left_mult_mps.at(b1, local);

            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1,b2);

        for (std::size_t op_index = 0; op_index < access.size(); ++op_index)
        {
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op(op_index);
            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);

            int a = mpo.left_spin(b1).get(), k = W.spin().get(), ap = mpo.right_spin(b2).get();

            for (size_t t_block = 0; t_block < T.n_blocks(); ++t_block) {

                charge lc = T.basis().left_charge(t_block);
                charge rc = T.basis().right_charge(t_block);

                const_iterator it = ket_basis.left_lower_bound(rc); //ket_basis comes transposed!
                charge mc = it->rc;

                for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                {
                    charge phys_in = W.basis().left_charge(w_block);
                    charge phys_out = W.basis().right_charge(w_block);

                    charge out_r_charge = SymmGroup::fuse(rc, phys_in);
                    size_t rb = right_i.position(out_r_charge);
                    if (rb == right_i.size()) continue;

                    charge out_l_charge = SymmGroup::fuse(lc, phys_out);
                    if (!::SU2::triangle(SymmGroup::spin(out_r_charge), ap, SymmGroup::spin(out_l_charge))) continue;
                    if (!right_i.has(out_l_charge)) continue; // can also probe out_left_i, but right_i has the same charges

                    size_t r_size = right_i[rb].second;

                    size_t o = ret.find_block(out_l_charge, out_r_charge);
                    if ( o == ret.n_blocks() )
                        o = ret.insert_block(Matrix(out_left_i.size_of_block(out_l_charge),r_size), out_l_charge, out_r_charge);

                    int i = SymmGroup::spin(lc), ip = SymmGroup::spin(out_l_charge);
                    int j = SymmGroup::spin(mc), jp = SymmGroup::spin(out_r_charge);
                    int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                    typename Matrix::value_type couplings[4];
                    ::SU2::set_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip, access.scale(op_index), couplings);

                    size_t in_right_offset = in_right_pb(phys_in, out_r_charge);
                    size_t out_left_offset = out_left_pb(phys_out, lc);
                    size_t l_size = T.basis().left_size(t_block);

                    detail::lbtm<Matrix, SymmGroup>(T[t_block], ret[o], W, in_right_offset, out_left_offset, l_size, r_size, w_block, couplings);
                } // wblock
            } // lblock
        } // op_index
        } // b1
    }

    template <class Matrix, class SymmGroup>
    struct task_capsule
    {
        typedef typename SymmGroup::charge charge;
        typedef typename Matrix::value_type value_type;
        typedef detail::micro_task<value_type> micro_task;
        typedef std::map<std::pair<charge, charge>, std::vector<micro_task>, compare_pair<std::pair<charge, charge> > > map_t;

        map_t tasks;        
    };

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void rbtm_tasks(size_t b1,
                    MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> const & right_mult_mps,
                    MPOTensor<Matrix, SymmGroup> const & mpo,
                    DualIndex<SymmGroup> const & ket_basis,
                    Index<SymmGroup> const & left_i,
                    Index<SymmGroup> const & out_right_i,
                    ProductBasis<SymmGroup> const & in_left_pb,
                    ProductBasis<SymmGroup> const & out_right_pb,
                    task_capsule<Matrix, SymmGroup> & tasks_cap)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename SymmGroup::charge charge;
        typedef typename Matrix::value_type value_type;

        typedef typename task_capsule<Matrix, SymmGroup>::micro_task micro_task;

        row_proxy row_b1 = mpo.row(b1);
        for (typename row_proxy::const_iterator row_it = row_b1.begin(); row_it != row_b1.end(); ++row_it) {
            index_type b2 = row_it.index();

            block_matrix<Matrix, SymmGroup> const & T = right_mult_mps.at(b2);

            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1,b2);

            for (size_t op_index = 0; op_index < access.size(); ++op_index)
            {
                typename operator_selector<Matrix, SymmGroup>::type const & W = access.op(op_index);
                int a = mpo.left_spin(b1).get(), k = W.spin().get(), ap = mpo.right_spin(b2).get();

                for (size_t t_block = 0; t_block < T.n_blocks(); ++t_block){

                    charge lc = T.basis().left_charge(t_block);
                    charge rc = T.basis().right_charge(t_block);

                    const_iterator it = ket_basis.left_lower_bound(lc);
                    charge mc = it->rc;

                    for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                    {   
                        charge phys_in = W.basis().left_charge(w_block);
                        charge phys_out = W.basis().right_charge(w_block);

                        charge out_l_charge = SymmGroup::fuse(lc, -phys_in);
                        size_t lb = left_i.position(out_l_charge);
                        if (lb == left_i.size()) continue;

                        charge out_r_charge = SymmGroup::fuse(rc, -phys_out);
                        if (!::SU2::triangle(SymmGroup::spin(out_l_charge), a, SymmGroup::spin(out_r_charge))) continue;
                        if (!left_i.has(out_r_charge)) continue;

                        size_t l_size = left_i[lb].second; 

                        std::vector<micro_task> & otasks = tasks_cap.tasks[std::make_pair(out_l_charge, out_r_charge)];

                        int i = SymmGroup::spin(out_r_charge), ip = SymmGroup::spin(rc);
                        int j = SymmGroup::spin(out_l_charge), jp = SymmGroup::spin(mc);
                        int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                        typename Matrix::value_type couplings[4];
                        ::SU2::set_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip, access.scale(op_index), couplings);

                        size_t in_left_offset = in_left_pb(phys_in, out_l_charge);
                        size_t out_right_offset = out_right_pb(phys_out, rc);
                        size_t r_size = T.basis().right_size(t_block);

                        micro_task tpl; tpl.l_size = l_size; tpl.stripe = num_rows(T[t_block]); tpl.b2 = b2; tpl.k = t_block;

                        unsigned short r_size_cache = 16384 / (l_size * W.basis().right_size(w_block)); // 128 KB
                        for (unsigned short slice = 0; slice < r_size/r_size_cache; ++slice)
                        {
                            size_t right_offset_cache = slice * r_size_cache;
                            detail::op_iterate<Matrix, SymmGroup>(W, w_block, couplings, otasks, tpl, 
                                                                  //&T[t_block](in_left_offset, right_offset_cache),
                                                                  in_left_offset + tpl.stripe * right_offset_cache,
                                                                  r_size_cache, r_size,
                                                                  out_right_offset + right_offset_cache);
                        }

                        unsigned short r_size_remain = r_size % r_size_cache;
                        unsigned short right_offset_remain = r_size - r_size_remain;
                        if (r_size_remain == 0) continue;

                        detail::op_iterate<Matrix, SymmGroup>(W, w_block, couplings, otasks, tpl,
                                                              //&T[t_block](in_left_offset, right_offset_remain),
                                                              in_left_offset + tpl.stripe * right_offset_remain,
                                                              r_size_remain, r_size,
                                                              out_right_offset + right_offset_remain);
                } // wblock
                } // ket block
            } // op_index
        } // b2
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void rbtm_axpy(task_capsule<Matrix, SymmGroup> & tasks_cap, block_matrix<Matrix, SymmGroup> & ret,
                   Index<SymmGroup> const & out_right_i,
                   MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> const & t)
    {
        typedef typename Matrix::value_type value_type;
        typedef typename task_capsule<Matrix, SymmGroup>::map_t map_t;
        typedef typename task_capsule<Matrix, SymmGroup>::micro_task micro_task;

        map_t & tasks = tasks_cap.tasks;
        for (typename map_t::iterator it = tasks.begin(); it != tasks.end(); ++it)
        {
            std::vector<micro_task> & otasks = it->second;
            std::sort(otasks.begin(), otasks.end(), detail::task_compare<value_type>()); 

            if (otasks.size() == 0) continue;
            Matrix buf(otasks[0].l_size, out_right_i.size_of_block(it->first.second));

            for (typename std::vector<micro_task>::const_iterator it2 = otasks.begin(); it2 != otasks.end(); ++it2)
                detail::task_axpy(*it2, &buf(0,0), &t.at(it2->b2)[it2->k](0,0) + it2->in_offset);

            ret.insert_block(buf, it->first.first, it->first.second);
        }
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void rbtm_kernel(size_t b1,
                     block_matrix<Matrix, SymmGroup> & ret,
                     Boundary<OtherMatrix, SymmGroup> const & left,
                     MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> const & right_mult_mps,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     DualIndex<SymmGroup> const & ket_basis,
                     Index<SymmGroup> const & left_i,
                     Index<SymmGroup> const & out_right_i,
                     ProductBasis<SymmGroup> const & in_left_pb,
                     ProductBasis<SymmGroup> const & out_right_pb)
    {
        task_capsule<Matrix, SymmGroup> tasks_cap;

        rbtm_tasks(b1, right_mult_mps, mpo, ket_basis, left_i, out_right_i, in_left_pb, out_right_pb, tasks_cap);
        rbtm_axpy(tasks_cap, ret, out_right_i, right_mult_mps);

        right_mult_mps.free(b1);

    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void charge_gemm(Matrix const & A, OtherMatrix const & B, block_matrix<OtherMatrix, SymmGroup> & C,
                     typename SymmGroup::charge rc, typename Matrix::value_type scale)
    {
        size_t c_block = C.find_block(rc, rc);
        if (c_block == C.n_blocks())
               c_block = C.insert_block(OtherMatrix(num_rows(A), num_cols(B)), rc, rc);

        boost::numeric::bindings::blas::gemm(scale, A, B, typename Matrix::value_type(1), C[c_block]); 
    }

    template<class Matrix, class OtherMatrix, class TVMatrix, class SymmGroup>
    void rbtm_axpy_gemm(size_t b1, task_capsule<Matrix, SymmGroup> & tasks_cap, block_matrix<Matrix, SymmGroup> & prod,
                        Index<SymmGroup> const & out_right_i, Boundary<OtherMatrix, SymmGroup> const & left,
                        MPOTensor<Matrix, SymmGroup> const & mpo,
                        block_matrix<TVMatrix, SymmGroup> const & left_b1,
                        MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> const & t)
    {
        typedef typename Matrix::value_type value_type;
        typedef typename task_capsule<Matrix, SymmGroup>::map_t map_t;
        typedef typename task_capsule<Matrix, SymmGroup>::micro_task micro_task;
        typedef typename SymmGroup::charge charge;

        std::vector<value_type> phases = (mpo.herm_info.left_skip(b1)) ? ::contraction::common::conjugate_phases(left_b1, mpo, b1, true, false) :
                                                                         std::vector<value_type>(left_b1.n_blocks(),1.);
        map_t & tasks = tasks_cap.tasks;
        for (typename map_t::iterator it = tasks.begin(); it != tasks.end(); ++it)
        {
            std::vector<micro_task> & otasks = it->second;
            std::sort(otasks.begin(), otasks.end(), detail::task_compare<value_type>()); 

            if (otasks.size() == 0) continue;
            Matrix buf(otasks[0].l_size, out_right_i.size_of_block(it->first.second));

            size_t k = left_b1.basis().position(it->first.second, it->first.first); if (k == left_b1.basis().size()) continue;

            for (typename std::vector<micro_task>::const_iterator it2 = otasks.begin(); it2 != otasks.end(); ++it2)
                detail::task_axpy(*it2, &buf(0,0), &t.at(it2->b2)[it2->k](0,0) + it2->in_offset);

            charge_gemm(left_b1[k], buf, prod, it->first.second, phases[k]);
        }
    }
} // namespace SU2
} // namespace contraction

#endif
