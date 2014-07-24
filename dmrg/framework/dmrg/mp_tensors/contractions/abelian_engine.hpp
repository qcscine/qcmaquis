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

#ifndef ABELIAN_ENGINE
#define ABELIAN_ENGINE

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/block_matrix/indexing.h"

#include "dmrg/mp_tensors/contractions/engine.h"
#include "dmrg/mp_tensors/contractions/abelian/functors.h"

// implementation includes see below

namespace contraction {

template <class Matrix, class OtherMatrix, class SymmGroup>
class AbelianEngine : public Engine<Matrix, OtherMatrix, SymmGroup>
{
public:
    AbelianEngine() {}

    virtual block_matrix<OtherMatrix, SymmGroup>
    overlap_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                      MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                      block_matrix<OtherMatrix, SymmGroup> const & left,
                      block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
    {
        return contraction::overlap_left_step<Matrix, OtherMatrix, SymmGroup, AbelianGemms>(bra_tensor, ket_tensor, left, localop);
    }


    virtual block_matrix<OtherMatrix, SymmGroup>
    overlap_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                       MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                       block_matrix<OtherMatrix, SymmGroup> const & right,
                       block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
    {
        return contraction::overlap_right_step<Matrix, OtherMatrix, SymmGroup, AbelianGemms>(bra_tensor, ket_tensor, right, localop);
    }

    virtual Boundary<Matrix, SymmGroup>
    left_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                             Boundary<OtherMatrix, SymmGroup> const & left,
                             MPOTensor<Matrix, SymmGroup> const & mpo,
                             Index<SymmGroup> const * in_low = NULL)
    {
        return contraction::left_boundary_tensor_mpo<Matrix, OtherMatrix, SymmGroup, AbelianGemms, lbtm_functor>
               (mps, left, mpo, in_low);
    }

    virtual Boundary<Matrix, SymmGroup>
    right_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                              Boundary<OtherMatrix, SymmGroup> const & right,
                              MPOTensor<Matrix, SymmGroup> const & mpo,
                              Index<SymmGroup> const * in_low = NULL)
    {
        return contraction::right_boundary_tensor_mpo<Matrix, OtherMatrix, SymmGroup, AbelianGemms, rbtm_functor>
               (mps, right, mpo, in_low);
    }

    virtual Boundary<OtherMatrix, SymmGroup>
    overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          Boundary<OtherMatrix, SymmGroup> const & left,
                          MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        return contraction::overlap_mpo_left_step<Matrix, OtherMatrix, SymmGroup, AbelianGemms, lbtm_functor>
               (bra_tensor, ket_tensor, left, mpo);
    }

    virtual Boundary<OtherMatrix, SymmGroup>
    overlap_mpo_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                           MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                           Boundary<OtherMatrix, SymmGroup> const & right,
                           MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        return contraction::overlap_mpo_right_step<Matrix, OtherMatrix, SymmGroup, AbelianGemms, rbtm_functor>
               (bra_tensor, ket_tensor, right, mpo);
    }

    virtual MPSTensor<Matrix, SymmGroup>
    site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
                Boundary<OtherMatrix, SymmGroup> const & left,
                Boundary<OtherMatrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo);

    virtual std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
    predict_new_state_l2r_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                MPOTensor<Matrix, SymmGroup> const & mpo,
                                Boundary<OtherMatrix, SymmGroup> const & left,
                                Boundary<OtherMatrix, SymmGroup> const & right,
                                double alpha, double cutoff, std::size_t Mmax)
    {
        return contraction::predict_new_state_l2r_sweep<Matrix, OtherMatrix, SymmGroup, AbelianGemms, lbtm_functor>
              (mps, mpo, left, right, alpha, cutoff, Mmax);
    }

    virtual MPSTensor<Matrix, SymmGroup>
    predict_lanczos_l2r_sweep(MPSTensor<Matrix, SymmGroup> B,
                              MPSTensor<Matrix, SymmGroup> const & psi,
                              MPSTensor<Matrix, SymmGroup> const & A)
    {
        return contraction::predict_lanczos_l2r_sweep<Matrix, SymmGroup, AbelianGemms>(B, psi, A);
    }

    virtual std::pair<MPSTensor<Matrix, SymmGroup>, truncation_results>
    predict_new_state_r2l_sweep(MPSTensor<Matrix, SymmGroup> const & mps,
                                MPOTensor<Matrix, SymmGroup> const & mpo,
                                Boundary<OtherMatrix, SymmGroup> const & left,
                                Boundary<OtherMatrix, SymmGroup> const & right,
                                double alpha, double cutoff, std::size_t Mmax)
    {
        return contraction::predict_new_state_r2l_sweep<Matrix, OtherMatrix, SymmGroup, AbelianGemms, rbtm_functor>
               (mps, mpo, left, right, alpha, cutoff, Mmax);             
    }

    virtual MPSTensor<Matrix, SymmGroup>
    predict_lanczos_r2l_sweep(MPSTensor<Matrix, SymmGroup> B,
                              MPSTensor<Matrix, SymmGroup> const & psi,
                              MPSTensor<Matrix, SymmGroup> const & A)
    {
        return contraction::predict_lanczos_r2l_sweep<Matrix, SymmGroup, AbelianGemms>(B, psi, A);
    }
};


template<class Matrix, class OtherMatrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup>
AbelianEngine<Matrix, OtherMatrix, SymmGroup>::
site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
            Boundary<OtherMatrix, SymmGroup> const & left,
            Boundary<OtherMatrix, SymmGroup> const & right,
            MPOTensor<Matrix, SymmGroup> const & mpo)
{
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

    std::vector<block_matrix<Matrix, SymmGroup> > t
        = boundary_times_mps<Matrix, OtherMatrix, SymmGroup, gemm_trim_left_functor>(ket_tensor, left, mpo);

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
        lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.data().basis(), right_i, out_left_i, in_right_pb, out_left_pb);
        contr_grid.multiply_column(b2, right[b2]);
    });
    t.clear();
    ambient::sync();

    swap(ret.data(), contr_grid.reduce());
    ambient::sync();
#else
    omp_for(index_type b2, range<index_type>(0,loop_max), {
        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);

        lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.data().basis(), right_i, out_left_i, in_right_pb, out_left_pb);
        block_matrix<Matrix, SymmGroup> tmp;
        gemm(contr_grid(0,0), right[b2], tmp);
        contr_grid(0,0).clear();
        omp_critical
        for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
            ret.data().match_and_add_block(tmp[k], tmp.left_basis_charge(k), tmp.right_basis_charge(k));
    });
#endif
    return ret;
}

} // namespace contractions

#endif
