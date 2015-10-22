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

    //template<class Matrix, class OtherMatrix, class SymmGroup>
    //void h_diag(size_t b2,
    //            ContractionGrid<typename alps::numeric::associated_diagonal_matrix<Matrix>::type, SymmGroup>& contr_grid,
    //            Boundary<OtherMatrix, SymmGroup> const & left,
    //            std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
    //            MPOTensor<Matrix, SymmGroup> const & mpo,
    //            DualIndex<SymmGroup> const & ket_basis,
    //            Index<SymmGroup> const & right_i,
    //            Index<SymmGroup> const & out_left_i,
    //            ProductBasis<SymmGroup> const & in_right_pb,
    //            ProductBasis<SymmGroup> const & out_left_pb)
    //{
    template<class Matrix, class OtherMatrix, class SymmGroup>
    typename Matrix::value_type
    h_diag(size_t b2,
           Boundary<OtherMatrix, SymmGroup> const & left,
           MPOTensor<Matrix, SymmGroup> const & mpo,
           DualIndex<SymmGroup> const & ket_basis,
           Index<SymmGroup> const & right_i,
           ProductBasis<SymmGroup> const & in_right_pb,
           ProductBasis<SymmGroup> const & out_left_pb,
           typename SymmGroup::charge lc,
           typename SymmGroup::charge rc,
           typename SymmGroup::charge pc,
           size_t in_phys_offset,
           size_t row_i)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename SymmGroup::charge charge;

        typename Matrix::value_type ret = 0.0;

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();

            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1,b2);

        for (std::size_t op_index = 0; op_index < access.size(); ++op_index)
        {
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op(op_index);
            int a = mpo.left_spin(b1).get(), k = W.spin().get(), ap = mpo.right_spin(b2).get();

            size_t l_block = left[b1].find_block(lc, lc);
            if (l_block == left[b1].n_blocks()) continue;

            //for (size_t lblock = 0; lblock < left[b1].n_blocks(); ++lblock)
            //{
            for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
            {
                charge phys_in = W.basis().left_charge(w_block);
                charge phys_out = W.basis().right_charge(w_block);
                if (phys_in != pc || phys_in != phys_out) continue;

                charge out_r_charge = SymmGroup::fuse(rc, phys_in);
                if (!right_i.has(out_r_charge)) continue;

                charge out_l_charge = SymmGroup::fuse(lc, phys_out);
                if (!right_i.has(out_l_charge)) continue;

                int i = SymmGroup::spin(lc), ip = SymmGroup::spin(out_l_charge);
                int j = SymmGroup::spin(lc), jp = SymmGroup::spin(out_r_charge);
                int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                typename Matrix::value_type prefactor = sqrt((ip+1.)*(j+1.)/((i+1.)*(jp+1.))) * access.scale(op_index);
                typename Matrix::value_type couplings[2];
                couplings[0] = prefactor * ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                couplings[1] = prefactor * ::SU2::mod_coupling(j, 2,     jp, a,k,ap, i, 2,      ip);

                size_t phys_s1 = W.basis().left_size(w_block);
                size_t phys_s2 = W.basis().right_size(w_block);

                typedef typename SparseOperator<Matrix, SymmGroup>::const_iterator block_iterator;
                std::pair<block_iterator, block_iterator> blocks = W.get_sparse().block(w_block);

                for( block_iterator it = blocks.first; it != blocks.second; ++it)
                {
                    std::size_t ss1 = it->row;
                    std::size_t ss2 = it->col;
                    if (in_phys_offset != ss1 || ss1 != ss2) continue;

                    std::size_t rspin = it->row_spin;
                    std::size_t cspin = it->col_spin;
                    std::size_t casenr = 0;
                    if (rspin == 2 && cspin == 2) casenr = 1;

                    typename Matrix::value_type alfa_t = it->coefficient * couplings[casenr];
                    ret += alfa_t * left[b1][l_block](row_i, row_i);
                }

                } // wblock
            //} // lblock
        } // op_index
        } // b1
        return ret;
    }

} // namespace SU2
} // namespace contraction

#endif
