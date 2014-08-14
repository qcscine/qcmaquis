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

#ifndef CONTRACTIONS_SITE_HAMIL_H
#define CONTRACTIONS_SITE_HAMIL_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/block_matrix/indexing.h"

namespace contraction {

    template<class Matrix, class OtherMatrix, class SymmGroup>
    MPSTensor<Matrix, SymmGroup>
    site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
                Boundary<OtherMatrix, SymmGroup> const & left,
                Boundary<OtherMatrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

        std::vector<block_matrix<Matrix, SymmGroup> > t = boundary_times_mps(ket_tensor, left, mpo);

        Index<SymmGroup> const & physical_i = ket_tensor.site_dim(),
                               & left_i = ket_tensor.row_dim(),
                               & right_i = ket_tensor.col_dim(),
                                 out_left_i = physical_i * left_i;
        ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> in_right_pb(physical_i, right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));

        MPSTensor<Matrix, SymmGroup> ret;
        ret.phys_i = ket_tensor.site_dim(); ret.left_i = ket_tensor.row_dim(); ret.right_i = ket_tensor.col_dim();
        index_type loop_max = mpo.col_dim();

#ifdef USE_AMBIENT
        {
            block_matrix<Matrix, SymmGroup> empty;
            swap(ket_tensor.data(), empty); // deallocating mpstensor before exiting the stack
        }
        ambient::sync(); 
        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, left.aux_dim(), mpo.col_dim());
        contr_grid.hint_left(t);
        contr_grid.hint_right(right);

        parallel_for(index_type b2, range<index_type>(0,loop_max), {
            lbtm_kernel(b2, contr_grid, left, t, mpo, physical_i, right_i, out_left_i, in_right_pb, out_left_pb);
            contr_grid.multiply_column(b2, right[b2]);
        });
        t.clear();
        ambient::sync();
         
        swap(ret.data(), contr_grid.reduce());
        ambient::sync();
#else
        omp_for(index_type b2, range<index_type>(0,loop_max), {
            ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);
            lbtm_kernel(b2, contr_grid, left, t, mpo, physical_i, right_i, out_left_i, in_right_pb, out_left_pb);
            block_matrix<Matrix, SymmGroup> tmp;
            gemm(contr_grid(0,0), right[b2], tmp);
            contr_grid(0,0).clear();
            omp_critical
            for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
                ret.data().match_and_add_block(tmp[k], tmp.basis().left_charge(k), tmp.basis().right_charge(k));
        });
#endif
        return ret;
    }
}
#endif
