/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Alexandr Kosenkov <akosenko@phys.ethz.ch>
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

#ifndef CONTRACTIONS_COMMON_BOUNDARY_TIMES_MPS_H
#define CONTRACTIONS_COMMON_BOUNDARY_TIMES_MPS_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/block_matrix/indexing.h"

namespace contraction {

    template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
    std::vector<block_matrix<OtherMatrix, SymmGroup> >
    boundary_times_mps(MPSTensor<Matrix, SymmGroup> const & mps,
                       Boundary<OtherMatrix, SymmGroup> const & left,
                       MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        std::vector<block_matrix<OtherMatrix, SymmGroup> > ret(left.aux_dim());
        int loop_max = left.aux_dim();
        {
            select_proc(storage::actor_t::common);
            mps.make_right_paired();
            storage::hint(mps);
        }
        parallel_for(int b1, range(0,loop_max), {
            select_proc(ambient::scope::permute(b1, mpo.placement_l));
            Gemm()(transpose(left[b1]), mps.data(), ret[b1]);
        });
        return ret;
    }

    template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
    std::vector<block_matrix<OtherMatrix, SymmGroup> >
    mps_times_boundary(MPSTensor<Matrix, SymmGroup> const & mps,
                       Boundary<OtherMatrix, SymmGroup> const & right,
                       MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        std::vector<block_matrix<OtherMatrix, SymmGroup> > ret(right.aux_dim());
        int loop_max = right.aux_dim();
        {
            select_proc(storage::actor_t::common);
            mps.make_left_paired();
            storage::hint(mps);
        }
        parallel_for(int b2, range(0,loop_max), {
            select_proc(ambient::scope::permute(b2, mpo.placement_r));
            Gemm()(mps.data(), right[b2], ret[b2]);
        });
        return ret;
    }

}
#endif
