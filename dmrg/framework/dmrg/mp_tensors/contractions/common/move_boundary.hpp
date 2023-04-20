/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef ENGINE_COMMON_MOVE_BOUNDARY_H
#define ENGINE_COMMON_MOVE_BOUNDARY_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/contractions/detail/memsave.hpp"

namespace contraction {
namespace common {

// output/input: left_i for bra_tensor, right_i for ket_tensor
template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
static block_matrix<OtherMatrix, SymmGroup>
overlap_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                  MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                  block_matrix<OtherMatrix, SymmGroup> const & left,
                  block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
{
    if (localop != NULL)
        throw std::runtime_error("Not implemented!");
    assert(ket_tensor.phys_i == bra_tensor.phys_i);
    block_matrix<OtherMatrix, SymmGroup> t1;
    block_matrix<Matrix, SymmGroup> t3;
    ket_tensor.make_right_paired();
    typename Gemm::gemm()(left, ket_tensor.data(), t1);
    reshape_right_to_left_new(ket_tensor.site_dim(), bra_tensor.row_dim(), ket_tensor.col_dim(),
                              t1, t3);
    bra_tensor.make_left_paired();
    typename Gemm::gemm()(transpose(conjugate(bra_tensor.data())), t3, t1);
    return t1;
    // original:
    // t3 = transpose(t3);
    // gemm(t3, t2, t1);
    // return transpose(t1);
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm>
static block_matrix<OtherMatrix, SymmGroup>
overlap_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                   MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                   block_matrix<OtherMatrix, SymmGroup> const & right,
                   block_matrix<OtherMatrix, SymmGroup> * localop = NULL)
{
    if (localop != NULL)
        throw std::runtime_error("Not implemented!");
    assert(ket_tensor.phys_i == bra_tensor.phys_i);
    ket_tensor.make_left_paired();
    block_matrix<OtherMatrix, SymmGroup> t1;
    block_matrix<Matrix, SymmGroup> t3;
    typename Gemm::gemm()(ket_tensor.data(), transpose(right), t1);
    reshape_left_to_right_new(ket_tensor.site_dim(), ket_tensor.row_dim(), bra_tensor.col_dim(), t1, t3);
    bra_tensor.make_right_paired();
    typename Gemm::gemm()(conjugate(bra_tensor.data()), transpose(t3), t1);
    return t1;
}

// note: this function changes the internal structure of Boundary,
//       each block is transposed
template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
static Boundary<Matrix, SymmGroup>
left_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                         Boundary<OtherMatrix, SymmGroup> const & left,
                         MPOTensor<Matrix, SymmGroup> const & mpo,
                         Index<SymmGroup> const * in_low = NULL)
{
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    if (in_low == NULL)
        in_low = &mps.row_dim();
    //std::vector<block_matrix<Matrix, SymmGroup> > t
    Index<SymmGroup> physical_i = mps.site_dim(), left_i = *in_low, right_i = mps.col_dim(),
                                  out_left_i = physical_i * left_i;
    BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemm> t(mps, left, mpo, left_i);
    ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
    ProductBasis<SymmGroup> in_right_pb(physical_i, right_i,
                            boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                    -boost::lambda::_1, boost::lambda::_2));
    index_type loop_max = mpo.col_dim();
    Boundary<Matrix, SymmGroup> ret;
    ret.resize(mpo.col_dim());
    omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);
        Kernel()(b2, contr_grid, left, t, mpo, mps.data().basis(), mps.data().basis(), 
                 right_i, out_left_i, in_right_pb, out_left_pb, true);
        swap(ret[b2], contr_grid(0,0));
    });
    return ret;
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
static Boundary<Matrix, SymmGroup>
right_boundary_tensor_mpo(MPSTensor<Matrix, SymmGroup> mps,
                          Boundary<OtherMatrix, SymmGroup> const & right,
                          MPOTensor<Matrix, SymmGroup> const & mpo,
                          Index<SymmGroup> const * in_low = NULL)
{
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    if (in_low == NULL)
        in_low = &mps.col_dim();
    contraction::common::MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, Gemm> t(mps, right, mpo);
    Index<SymmGroup> physical_i = mps.site_dim(), left_i = mps.row_dim(), right_i = *in_low,
                     out_right_i = adjoin(physical_i) * right_i;
    ProductBasis<SymmGroup> in_left_pb(physical_i, left_i);
    ProductBasis<SymmGroup> out_right_pb(physical_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
    Boundary<Matrix, SymmGroup> ret;
    ret.resize(mpo.row_dim());
    index_type loop_max = mpo.row_dim();
    omp_for(index_type b1, parallel::range<index_type>(0,loop_max), {
        // parallel::guard group(scheduler(b1), parallel::groups_granularity);
        Kernel()(b1, ret[b1], right, t, mpo, mps.data().basis(), mps.data().basis(),
                 left_i, out_right_i, in_left_pb, out_right_pb, true);
    });
    return ret;
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
static Boundary<OtherMatrix, SymmGroup>
overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor, MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                      Boundary<OtherMatrix, SymmGroup> const & left, MPOTensor<Matrix, SymmGroup> const & mpo,
                      bool isHermitian=true)
{
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    bra_tensor.make_right_paired();
    Index<SymmGroup> braBasis = bra_tensor.data().left_basis();
    MPSTensor<Matrix, SymmGroup> ket_cpy = ket_tensor;
    BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemm> t(ket_cpy, left, mpo, braBasis, isHermitian);
    Index<SymmGroup> const & left_i = bra_tensor.row_dim();
    Index<SymmGroup> right_i = ket_tensor.col_dim();
    Index<SymmGroup> bra_right_i = bra_tensor.col_dim();
    Index<SymmGroup> out_left_i = bra_tensor.site_dim() * left_i;
    common_subset(out_left_i, bra_right_i);
    ProductBasis<SymmGroup> out_left_pb(bra_tensor.site_dim(), left_i);
    ProductBasis<SymmGroup> in_right_pb(ket_tensor.site_dim(), right_i,
                            boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                    -boost::lambda::_1, boost::lambda::_2));
    index_type loop_max = mpo.col_dim();
    DualIndex<SymmGroup> bra_basis = bra_tensor.data().basis();
    bra_tensor.make_left_paired();
    block_matrix<Matrix, SymmGroup> bra_conj = conjugate(bra_tensor.data());
    DualIndex<SymmGroup> ket_basis_transpose = ket_cpy.data().basis();
    for (std::size_t i = 0; i < ket_basis_transpose.size(); ++i) {
        std::swap(ket_basis_transpose[i].lc, ket_basis_transpose[i].rc);
        std::swap(ket_basis_transpose[i].ls, ket_basis_transpose[i].rs);
    }
    Boundary<Matrix, SymmGroup> ret;
    ret.resize(loop_max);
    omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
        if (mpo.herm_info.right_skip(b2) && isHermitian)
            continue;
        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);
        Kernel()(b2, contr_grid, left, t, mpo, ket_basis_transpose, bra_basis, right_i, out_left_i, in_right_pb, out_left_pb, isHermitian);
        typename Gemm::gemm()(transpose(contr_grid(0,0)), bra_conj, ret[b2], MPOTensor_detail::get_spin(mpo, b2, false));
    });
    /*
    // hermiticity check
    omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
      if (mpo.herm_info.right_skip(b2))
      {
        block_matrix<Matrix, SymmGroup> conj0 = ret[mpo.herm_info.right_conj(b2)];
        recover_conjugate(conj0, mpo, b2, false, false);
        block_matrix<Matrix, SymmGroup> conj = transpose(conj0);
        block_matrix<Matrix, SymmGroup> diff = conj-ret[b2];
        if ( diff.norm() > 1e-6)
        {
            maquis::cout << "diff: " << diff << std::endl;
            maquis::cout << "should:\n" << ret[b2];
            maquis::cout << "is:\n" << conj;
            exit(1);
        }
      }
    });
    */
    return ret;
}

template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
static Boundary<OtherMatrix, SymmGroup>
overlap_mpo_right_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor, MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                       Boundary<OtherMatrix, SymmGroup> const & right, MPOTensor<Matrix, SymmGroup> const & mpo,
                       bool isHermitian=true)
{
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    Index<SymmGroup> const & physical_i = ket_tensor.site_dim(),
                             right_i = bra_tensor.col_dim();
    MPSTensor<Matrix, SymmGroup> ket_cpy = ket_tensor;
    Index<SymmGroup> left_i = ket_tensor.row_dim(),
                     out_right_i = adjoin(physical_i) * right_i,
                     bra_left_i = bra_tensor.row_dim();
    bra_tensor.make_left_paired();
    Index<SymmGroup> indexForTrim = bra_tensor.data().right_basis();
    contraction::common::MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, Gemm> t(ket_cpy, right, mpo, indexForTrim, isHermitian);
    common_subset(out_right_i, bra_left_i);
    ProductBasis<SymmGroup> in_left_pb(physical_i, left_i);
    ProductBasis<SymmGroup> out_right_pb(physical_i, right_i,
                                         boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
    Boundary<Matrix, SymmGroup> ret;
    ret.resize(mpo.row_dim());
    //ket_tensor.make_right_paired();
    index_type loop_max = mpo.row_dim();
    bra_tensor.make_left_paired();
    DualIndex<SymmGroup> bra_basis = bra_tensor.data().basis();
    bra_tensor.make_right_paired();
    block_matrix<Matrix, SymmGroup> bra_conj = conjugate(bra_tensor.data());
    omp_for(index_type b1, parallel::range<index_type>(0,loop_max), {
        if (mpo.herm_info.left_skip(b1) && isHermitian)
            continue;
        Kernel()(b1, ret[b1], right, t, mpo, ket_cpy.data().basis(), bra_basis, left_i, out_right_i, in_left_pb, out_right_pb, isHermitian);
        block_matrix<Matrix, SymmGroup> tmp;
        typename Gemm::gemm()(ret[b1], transpose(bra_conj), tmp, MPOTensor_detail::get_spin(mpo, b1, true));
        //gemm(ret[b1], transpose(bra_conj), tmp, parallel::scheduler_size_indexed(ret[b1]));
        swap(ret[b1], tmp);
    });
    return ret;
}

// +-----------------------+
//  GENERATE_LEFT_MPO_BASIS
// +-----------------------+
template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
static Boundary<OtherMatrix, SymmGroup>
generate_left_mpo_basis(MPSTensor<Matrix, SymmGroup> const & bra_tensor,   // Bra MPSTensor object
                        MPSTensor<Matrix, SymmGroup> const & ket_tensor,   // Ket MPSTensor object
                        Boundary<OtherMatrix, SymmGroup> const & left,     // Left boundary
                        MPOTensor<Matrix, SymmGroup> const & mpo)          // MPOTensor object)
{
    // Types definition
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    MPSTensor<Matrix, SymmGroup> ket_cpy = ket_tensor;
    // Contracts the boundary with the MPS. Returns a vector of block_matrix objects.
    contraction::common::BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, Gemm> t(ket_cpy, left, mpo);
    // Here index are generated on the basis of the symmetry.
    // Consider that inside boundary_times_mps ket_cpy is made right paired, thus also the output
    // Will be right paired. However, the multiplication with the bra will contract both the physical
    // and the auxiliary index, thus the output has to be left paired.
    Index<SymmGroup> const & left_i = bra_tensor.row_dim();
    Index<SymmGroup> right_i = ket_tensor.col_dim();
    Index<SymmGroup> out_left_i = ket_tensor.site_dim() * left_i;
    common_subset(out_left_i, right_i);
    // The input product basis (in_right_pb) is computed as product of the site_dim of the ket_tensor object
    // by the its col_dim. In fact, the tensor in input is right_paired. The output is instead left_paired
    // (and the fusion has to be done with the minus sign).
    ProductBasis<SymmGroup> out_left_pb(ket_tensor.site_dim(), left_i);
    ProductBasis<SymmGroup> in_right_pb(ket_tensor.site_dim(), right_i,
                                        boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                            -boost::lambda::_1, boost::lambda::_2));
    index_type loop_max = mpo.col_dim();
    DualIndex<SymmGroup> ket_basis_transpose = ket_cpy.data().basis();
    for (std::size_t i = 0; i < ket_basis_transpose.size(); ++i) {
        std::swap(ket_basis_transpose[i].lc, ket_basis_transpose[i].rc);
        std::swap(ket_basis_transpose[i].ls, ket_basis_transpose[i].rs);
    }
    Boundary<Matrix, SymmGroup> ret;
    ret.resize(loop_max);
    // During the for cycle, which is run in parallel, a vector of block_matrix object is populated, which
    // is the grid attribute of the ContractionGrid object. The grid vector is then contracted back with
    // the bra tensor.
    omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
        if (mpo.herm_info.right_skip(b2))
            continue;
        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);
        Kernel()(b2, contr_grid, left, t, mpo, ket_basis_transpose, ket_basis_transpose, right_i, out_left_i,
                 in_right_pb, out_left_pb, true);
        // Final contraction with the MPS
        ret[b2] = contr_grid(0,0);
    });
    return ret;
};

// +------------------------+
//  GENERATE_RIGHT_MPO_BASIS
// +------------------------+
// Parallel of [GENERATE_LEFT_MPO_BASIS].
template<class Matrix, class OtherMatrix, class SymmGroup, class Gemm, class Kernel>
static Boundary<OtherMatrix, SymmGroup>
generate_right_mpo_basis(MPSTensor<Matrix, SymmGroup> const & bra_tensor, MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                         Boundary<OtherMatrix, SymmGroup> const & right, MPOTensor<Matrix, SymmGroup> const & mpo)
{
    // Types definition
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    // Contracts with the right boundary
    MPSTensor<Matrix, SymmGroup> ket_cpy = ket_tensor;
    contraction::common::MPSBoundaryProduct<Matrix, OtherMatrix, SymmGroup, Gemm> t(ket_cpy, right, mpo);
    // Generates the indexes
    Index<SymmGroup> const & physical_i = ket_tensor.site_dim();
    Index<SymmGroup> right_i = bra_tensor.col_dim();
    Index<SymmGroup> left_i = ket_tensor.row_dim();
    Index<SymmGroup> out_right_i = adjoin(physical_i) * right_i;
    common_subset(out_right_i, left_i);
    ProductBasis<SymmGroup> in_left_pb(physical_i, left_i);
    ProductBasis<SymmGroup> out_right_pb(physical_i, right_i, boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                             -boost::lambda::_1, boost::lambda::_2));
    // Prepares output
    Boundary<Matrix, SymmGroup> ret;
    ret.resize(mpo.row_dim());
    index_type loop_max = mpo.row_dim();
    // Main loop
    auto now = std::chrono::high_resolution_clock::now();
    omp_for(index_type b1, parallel::range<index_type>(0,loop_max), {
        if (mpo.herm_info.left_skip(b1))
            continue;
        Kernel()(b1, ret[b1], right, t, mpo, ket_cpy.data().basis(), ket_cpy.data().basis(), left_i, out_right_i, in_left_pb, out_right_pb, true);
    });
    return ret;
}

} // namespace common
} // namespace contraction

#endif
