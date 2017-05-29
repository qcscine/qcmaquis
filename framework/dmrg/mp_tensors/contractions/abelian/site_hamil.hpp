/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *                    Laboratory for Physical Chemistry, ETH Zurich
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

#ifndef ABELIAN_SITE_HAMIL
#define ABELIAN_SITE_HAMIL

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"

namespace contraction {

    template<class Matrix, class OtherMatrix, class SymmGroup, class SymmType>
    MPSTensor<Matrix, SymmGroup>
    Engine<Matrix, OtherMatrix, SymmGroup, SymmType>::
    site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
                Boundary<OtherMatrix, SymmGroup> const & left,
                Boundary<OtherMatrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

        std::vector<block_matrix<Matrix, SymmGroup> > t
            = common::boundary_times_mps<Matrix, OtherMatrix, SymmGroup, Gemms>(ket_tensor, left, mpo);
	    // Here indexes a grouped together
        Index<SymmGroup> const & physical_i = ket_tensor.site_dim(),
                               & left_i = ket_tensor.row_dim();
        Index<SymmGroup> right_i = ket_tensor.col_dim(),
                         out_left_i = physical_i * left_i;
        common_subset(out_left_i, right_i);
        ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> in_right_pb(physical_i, right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));
	    // Set the new MPSTensor
        MPSTensor<Matrix, SymmGroup> ret;
        ret.phys_i = ket_tensor.site_dim();
	    ret.left_i = ket_tensor.row_dim();
	    ret.right_i = ket_tensor.col_dim();
        index_type loop_max = mpo.col_dim();
        omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
            ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);
            abelian::lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.data().basis(), right_i, out_left_i, in_right_pb, out_left_pb);
            block_matrix<Matrix, SymmGroup> tmp;
            if (mpo.herm_info.right_skip(b2))
                gemm(contr_grid(0,0), transpose(right[mpo.herm_info.right_conj(b2)]), tmp);
            else
                gemm(contr_grid(0,0), right[b2], tmp);
            contr_grid(0,0).clear();
            parallel_critical
            for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
                ret.data().match_and_add_block(tmp[k], tmp.basis().left_charge(k), tmp.basis().right_charge(k));
        });
        return ret;
    }
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

            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op(0);
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

                            typedef typename SparseOperator<Matrix, SymmGroup>::const_iterator block_iterator;
                            std::pair<block_iterator, block_iterator> blocks = W.get_sparse().block(w_block);

                            for (block_iterator it = blocks.first; it != blocks.second; ++it)
                            {
                                std::size_t ss1 = it->row;
                                if (ss1 != it->col) continue;

                                for (size_t col_i = 0; col_i < right_i[block].second; ++col_i)
                                    std::transform(left_diagonal.begin(), left_diagonal.end(),
                                                   &ret[o](left_offset + ss1 * left_i[l].second, col_i),
                                                   &ret[o](left_offset + ss1 * left_i[l].second, col_i),
                                                   boost::lambda::_2 += boost::lambda::_1 );
                            }
                        } // wblock
                    } // phys_i s
                } // ket block
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
        // physical_i is the number of LOCAL basis functions
        // We left-pair the indexes.
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
} // namespace contraction

#endif
