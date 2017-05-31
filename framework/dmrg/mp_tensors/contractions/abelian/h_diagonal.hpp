/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*                    Laboratory for Physical Chemistry, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
*               2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef ABELIAN_H_DIAGONAL
#define ABELIAN_H_DIAGONAL

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"

namespace contraction {
//
// Function template for the calculation of the T matrix by assuming a diagonal MPO (see JCP 2015) for additional
// details. The quantity which is computed is:
//
//  -----  -----   -----          o_l o_l'      |    b_{l-1}          |\/| o_l'
//  \      \       \       \    /               |                     |  |
//  /      /       /        \/\/  b_{l-1} b_l   |___ a_{l-1} a_{l-1}' |  | a_{l-1}' a_l'
//  -----  -----   -----
//   o_l' a_{l-1}' b_{l-1}  |                                        |
//                          +----------------------------------------+
//                                   |
//                                   +--------------->    This is T (summed over b). In principle has 5
//                                                        indexes, but only 3 in this case due to symmetry.
//
// where W is the matrix involved in the definition of the MPO, L is the left boundary and M is the matrix
// involved in the definition of the MPS.
// We have < a_{l-1} a_l o_l | H | o_l' a_{l-1}' a_l' > . Taking the diagonal part means that we are
// taking o_l = o_l' , a_l = a_l' and a_{l-1} = a_{l-1}'. The only sum surviving is the one over b_{l-1}
// NOTE: the index b_l (column) is given as input and is fixed throughout the routine
//
//
    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    lbtm_diag_kernel(size_t b2,
                     Boundary<OtherMatrix, SymmGroup> const &left,
                     MPOTensor <Matrix, SymmGroup> const &mpo,
                     Index<SymmGroup> const &out_left_i,
                     Index<SymmGroup> const &left_i,
                     Index<SymmGroup> const &right_i,
                     Index<SymmGroup> const &phys_i,
                     ProductBasis<SymmGroup> const &left_pb) {
        //
        // Initialization
        // --------------
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename SymmGroup::charge charge;
        typedef typename Matrix::value_type value_type;
        block_matrix<Matrix, SymmGroup> ret;
        // Loop over the columns of the MPO (which is the dimension of the data of ref, i.e. the output matrix).
        // We sum over this dimension. The symmetry has no impact on the b indexes
        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();
            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1, b2);
            typename operator_selector<Matrix, SymmGroup>::type const &W = access.op(0);
            // Loop over a_l' (right dimension of ret). Pay attention to the symmetry
            for (size_t block = 0; block < right_i.size(); ++block) {
                charge in_charge = right_i[block].first;
                size_t o = ret.find_block(in_charge, in_charge);
                // Insert block with the same symmetry of the right index
                if (o == ret.n_blocks())
                    o = ret.insert_block(Matrix(out_left_i[block].second, right_i[block].second), in_charge, in_charge);
                // Loop over the physical indexes
                for (size_t s = 0; s < phys_i.size(); ++s) {
                    // Symmetry requirement : o_l*a_{l-1} == a_l'. In this way we determine a_{l-1}, which is also
                    // the index of the left boundary to be extracted.
                    charge phys_charge = phys_i[s].first;
                    size_t l = left_i.position(SymmGroup::fuse(in_charge, -phys_charge));
                    if (l == left_i.size()) continue;
                    charge lc = left_i[l].first;
                    size_t l_block = left[b1].find_block(lc, lc);
                    if (l_block == left[b1].n_blocks()) continue;
                    // copy the diagonal elements of the boundary into a vector
                    std::vector<value_type> left_diagonal(left_i[l].second);
                    // Diagonal is a method of the alps matrix class that gives two iterators, for the
                    // beginning and the end of the diagonal.
                    std::copy(left[b1][l_block].diagonal().first, left[b1][l_block].diagonal().second,
                              left_diagonal.begin());
                    // The overall size of (o_l*a_{l-1}) is given by the product
                    size_t left_offset = left_pb(phys_charge, lc);
                    // Loop over the basis of the MPO. Both the basis have to have the
                    // qn phys_charge
                    for (size_t w_block = 0; w_block < W.basis().size(); ++w_block) {
                        charge phys_in = W.basis().left_charge(w_block);
                        charge phys_out = W.basis().right_charge(w_block);
                        if (phys_charge != phys_in || phys_in != phys_out) continue;
                        typedef typename SparseOperator<Matrix, SymmGroup>::const_iterator block_iterator;
                        std::pair<block_iterator, block_iterator> blocks = W.get_sparse().block(w_block);
                        for (block_iterator it = blocks.first; it != blocks.second; ++it) {
                            std::size_t ss1 = it->row;
                            if (ss1 != it->col) continue;
                            for (size_t col_i = 0; col_i < right_i[block].second; ++col_i)
                                std::transform(left_diagonal.begin(), left_diagonal.end(),
                                               &ret[o](left_offset + ss1 * left_i[l].second, col_i),
                                               &ret[o](left_offset + ss1 * left_i[l].second, col_i),
                                               boost::lambda::_2 += boost::lambda::_1);
                        }
                    }
                }
            }
        }
        return ret;
    }

//
// Function template to contract T with the right boundaries.
//
//  -----     -----  o_l o_l      +--+  b_l       -+
//  \           |                 +--+             |-> this is has three index (one of which is sigma)
//  /           |    b_l a_{l-1}  |  \ a_l a_l    -+   as any other MPSTensor.
//  -----
//   a_l'
//
//
// where W is the matrix involved in the definition of the MPO, L is the left boundary and M is the matrix
// involved in the definition of the MPS.
// We have < a_{l-1} a_l o_l | H | o_l' a_{l-1}' a_l' > . Taking the diagonal part means that we are
// taking o_l = o_l' , a_l = a_l' and a_{l-1} = a_{l-1}'. The only sum surviving is the one over b_{l-1}
// NOTE: the index b_l (column) is given as input and is fixed throughout the routine
//
//
    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<Matrix, SymmGroup>
    diagonal_hamiltonian(Boundary<OtherMatrix, SymmGroup> const &left,
                         Boundary<OtherMatrix, SymmGroup> const &right,
                         MPOTensor <Matrix, SymmGroup> const &mpo,
                         MPSTensor<Matrix, SymmGroup> const &x) {
        //
        //
        // Initialization
        // --------------
        // Note that the phys index and the row are grouped together
        typedef typename SymmGroup::charge charge;
        Index<SymmGroup> const &physical_i = x.site_dim();
        Index<SymmGroup> right_i = x.col_dim();
        Index<SymmGroup> out_left_i = physical_i * x.row_dim();
        common_subset(out_left_i, right_i);
        ProductBasis<SymmGroup> out_left_pb(physical_i, x.row_dim());
        block_matrix<Matrix, SymmGroup> ret;
        // Loop over the b2 index, which is the "higher" index of the T tensor
        for (size_t b2 = 0; b2 < right.aux_dim(); ++b2) {
            block_matrix<Matrix, SymmGroup> lb2 = lbtm_diag_kernel(b2, left, mpo, out_left_i,
                                                                   x.row_dim(), x.col_dim(), physical_i,
                                                                   out_left_pb);
            for (size_t block = 0; block < lb2.n_blocks(); ++block) {
                charge in_r_charge = lb2.basis().right_charge(block);
                size_t rblock = right[b2].find_block(in_r_charge, in_r_charge);
                if (rblock != right[b2].n_blocks()) {
                    // Final contraction by
                    for (size_t c = 0; c < num_cols(lb2[block]); ++c)
                        std::transform(lb2[block].col(c).first, lb2[block].col(c).second, lb2[block].col(c).first,
                                       boost::lambda::_1 * right[b2][rblock](c, c));
                    ret.match_and_add_block(lb2[block], in_r_charge, in_r_charge);
                }
            }
        }
        return ret;
    }
}

#endif