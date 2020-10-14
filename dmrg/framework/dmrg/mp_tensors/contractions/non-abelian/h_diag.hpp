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

#ifndef CONTRACTIONS_SU2_H_DIAG_HPP
#define CONTRACTIONS_SU2_H_DIAG_HPP

#include "dmrg/block_matrix/symmetry/gsl_coupling.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

namespace contraction {
namespace SU2 {

    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    lbtm_diag_kernel(size_t b2,
                     Boundary<OtherMatrix, SymmGroup> const & left,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     Index<SymmGroup> const & out_left_i,
                     Index<SymmGroup> const & left_i,
                     Index<SymmGroup> const & right_i,
                     Index<SymmGroup> const & phys_i,
                     ProductBasis<SymmGroup> const & left_pb)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename SymmGroup::charge charge;
        typedef typename Matrix::value_type value_type;

        block_matrix<Matrix, SymmGroup> ret;

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();

            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1,b2);

        for (std::size_t op_index = 0; op_index < access.size(); ++op_index)
        {
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op(op_index);
            int a = mpo.left_spin(b1).get(), k = W.spin().get(), ap = mpo.right_spin(b2).get();

            // out_left_i and right_i have identical charges, but different sector sizes
            for (size_t block = 0; block < right_i.size(); ++block)
            {
                charge in_charge = right_i[block].first;

                size_t o = ret.find_block(in_charge, in_charge);
                if ( o == ret.n_blocks() )
                    o = ret.insert_block(Matrix(out_left_i[block].second, right_i[block].second), in_charge, in_charge);

                for (size_t s = 0; s < phys_i.size(); ++s)
                {
                    charge phys_charge = phys_i[s].first;
                    size_t l = left_i.position(SymmGroup::fuse(in_charge, -phys_charge));
                    if(l == left_i.size()) continue;
                    charge lc = left_i[l].first;
     
                    size_t l_block = left[b1].find_block(lc, lc);
                    if (l_block == left[b1].n_blocks()) continue;

                    // copy the diagonal elements of the boundary into a vector
                    std::vector<value_type> left_diagonal(left_i[l].second);
                    std::copy(left[b1][l_block].diagonal().first, left[b1][l_block].diagonal().second, left_diagonal.begin()); 

                    size_t left_offset = left_pb(phys_charge, lc);

                    for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                    {
                        charge phys_in = W.basis().left_charge(w_block);
                        charge phys_out = W.basis().right_charge(w_block);
                        if (phys_charge != phys_in || phys_in != phys_out) continue;

                        int i = SymmGroup::spin(lc), ip = SymmGroup::spin(in_charge);
                        int j = SymmGroup::spin(lc), jp = SymmGroup::spin(in_charge);
                        int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                        value_type prefactor = value_type(sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.)))) * access.scale(op_index);
                        value_type couplings[2];
                        couplings[0] = prefactor * (value_type)::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                        couplings[1] = prefactor * (value_type)::SU2::mod_coupling(j, 2,     jp, a,k,ap, i, 2,      ip);

                        typedef typename SparseOperator<Matrix, SymmGroup>::const_iterator block_iterator;
                        std::pair<block_iterator, block_iterator> blocks = W.get_sparse().block(w_block);

                        for (block_iterator it = blocks.first; it != blocks.second; ++it)
                        {
                            std::size_t ss1 = it->row;
                            if (ss1 != it->col) continue;

                            std::size_t rspin = it->row_spin;
                            std::size_t casenr = 0;
                            if (rspin == 2) casenr = 1;

                            typename Matrix::value_type alfa_t = it->coefficient * couplings[casenr];

                            // copy the diagonal multplied by alfa_t into the output result: ret(·, col_i) += alfa_t * diag(L)
                            for (size_t col_i = 0; col_i < right_i[block].second; ++col_i)
                                std::transform(left_diagonal.begin(), left_diagonal.end(),
                                               &ret[o](left_offset + ss1 * left_i[l].second, col_i),
                                               &ret[o](left_offset + ss1 * left_i[l].second, col_i),
                                               boost::lambda::_2 += boost::lambda::_1 * alfa_t);
                        }
                    } // wblock
                } // phys_i s
            } // ket block
        } // op_index
        } // b1
        return ret;
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    diagonal_hamiltonian(Boundary<OtherMatrix, SymmGroup> const & left,
                         Boundary<OtherMatrix, SymmGroup> const & right,
                         MPOTensor<Matrix, SymmGroup> const & mpo,                         
                         MPSTensor<Matrix, SymmGroup> const & x)
    {
        typedef typename SymmGroup::charge charge;

        Index<SymmGroup> const & physical_i = x.site_dim();
        Index<SymmGroup> right_i = x.col_dim(),
                         out_left_i = physical_i * x.row_dim();

        common_subset(out_left_i, right_i);
        ProductBasis<SymmGroup> out_left_pb(physical_i, x.row_dim());

        block_matrix<Matrix, SymmGroup> ret;
        for (size_t b2 = 0; b2 < right.aux_dim(); ++b2)
        {
            block_matrix<Matrix, SymmGroup> lb2 = lbtm_diag_kernel(b2, left, mpo, out_left_i,
                                                                   x.row_dim(), x.col_dim(), physical_i,
                                                                   out_left_pb);

            for (size_t block = 0; block < lb2.n_blocks(); ++block)
            {
                charge in_r_charge = lb2.basis().right_charge(block);
                size_t rblock = right[b2].find_block(in_r_charge, in_r_charge);
                if (rblock != right[b2].n_blocks())
                {
                    for (size_t c = 0; c < num_cols(lb2[block]); ++c)
                        std::transform(lb2[block].col(c).first, lb2[block].col(c).second, lb2[block].col(c).first,
                                       boost::lambda::_1 * right[b2][rblock](c,c));
                    ret.match_and_add_block(lb2[block], in_r_charge, in_r_charge);
                }
            }
        }
        return ret;
    }

} // namespace SU2
} // namespace contraction

#endif
