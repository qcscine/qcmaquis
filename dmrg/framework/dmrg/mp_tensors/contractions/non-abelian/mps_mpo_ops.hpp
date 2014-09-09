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

#ifndef CONTRACTIONS_SU2_MPS_MPO_OPS_HPP
#define CONTRACTIONS_SU2_MPS_MPO_OPS_HPP

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"

namespace SU2 {

    template<class Matrix, class SymmGroup>
    typename MPS<Matrix, SymmGroup>::scalar_type norm(MPS<Matrix, SymmGroup> const & mps)
    {
        std::size_t L = mps.length();

        block_matrix<Matrix, SymmGroup> left;
        left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);

        for(size_t i = 0; i < L; ++i) {
            select_proc(ambient::scope::balance(i,L));
            MPSTensor<Matrix, SymmGroup> cpy = mps[i];
            left = contraction::overlap_left_step<Matrix, Matrix, SymmGroup, SU2Gemms>(mps[i], cpy, left); // serial
        }

        return trace(left);
    }

    template<class Matrix, class SymmGroup>
    typename MPS<Matrix, SymmGroup>::scalar_type norm_r(MPS<Matrix, SymmGroup> const & mps)
    {
        std::size_t L = mps.length();

        block_matrix<Matrix, SymmGroup> right;
        right.insert_block(Matrix(1, 1, 1), mps[L-1].row_dim()[0].first, mps[L-1].row_dim()[0].first);

        for(int i = L-1; i >= 0 ; --i) {
            select_proc(ambient::scope::balance(i,L));
            MPSTensor<Matrix, SymmGroup> cpy = mps[i];
            right = contraction::overlap_right_step<Matrix, Matrix, SymmGroup, SU2Gemms>(mps[i], cpy, right); // serial
        }

        return trace(right);
    }

    template<class Matrix, class SymmGroup>
    double expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo,
                  boost::shared_ptr<contraction::Engine<Matrix, Matrix, SymmGroup> > contr,
                  int p1, int p2, std::vector<int> config)
    {
        assert(mpo.length() == mps.length());
        std::size_t L = mps.length();
        Boundary<Matrix, SymmGroup> left = mps.left_boundary();
        block_matrix<Matrix, SymmGroup> left_bm = mps.left_boundary()[0];

        for(size_t i = 0; i < L; ++i) {
            MPSTensor<Matrix, SymmGroup> cpy = mps[i];
            //left_bm = contraction::SU2::apply_operator(cpy, mps[i], left_bm, mpo[i], config, debug);
            left = contr->overlap_mpo_left_step(cpy, mps[i], left, mpo[i]);
        }

        return maquis::real(left[0].trace());
    }

    template<class Matrix, class SymmGroup>
    double expval_r(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo,
                    boost::shared_ptr<contraction::Engine<Matrix, Matrix, SymmGroup> > contr,
                    int p1, int p2, std::vector<int> config)
    {
        assert(mpo.length() == mps.length());
        std::size_t L = mps.length();
        Boundary<Matrix, SymmGroup> right = mps.right_boundary();

        for(int i = L-1; i >= 0; --i) {
            MPSTensor<Matrix, SymmGroup> cpy = mps[i];
            right = contr->overlap_mpo_right_step(cpy, mps[i], right, mpo[i]);
        }

        return maquis::real(right[0].trace());
    }

} // namespace SU2

#endif
