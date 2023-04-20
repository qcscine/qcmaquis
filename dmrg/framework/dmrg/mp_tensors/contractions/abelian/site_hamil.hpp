/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef ABELIAN_SITE_HAMIL
#define ABELIAN_SITE_HAMIL

namespace contraction {

template<class Matrix, class OtherMatrix, class SymmGroup, class SymmType>
MPSTensor<Matrix, SymmGroup> Engine<Matrix, OtherMatrix, SymmGroup, SymmType>::
  site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor,
              Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
              MPOTensor<Matrix, SymmGroup> const & mpo, bool isHermitian)
{
    return site_hamil2(ket_tensor, ket_tensor, left, right, mpo, isHermitian);
};


template<class Matrix, class OtherMatrix, class SymmGroup, class SymmType>
MPSTensor<Matrix, SymmGroup>
Engine<Matrix, OtherMatrix, SymmGroup, SymmType>::
site_hamil2(MPSTensor<Matrix, SymmGroup> ket_tensor, MPSTensor<Matrix, SymmGroup> const & bra_tensor,
            Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
            MPOTensor<Matrix, SymmGroup> const & mpo, bool isHermitian)
{
    // Type declaration
    using charge = typename SymmGroup::charge;
    using index_type = std::size_t;
    //
    Index<SymmGroup> const & physical_i = ket_tensor.site_dim();
    Index<SymmGroup> const & left_i = bra_tensor.row_dim();
    bra_tensor.make_right_paired();
    Index<SymmGroup> indexForTrim = bra_tensor.data().left_basis();
    contraction::common::BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, abelian::Gemms> t(ket_tensor, left, mpo, indexForTrim, isHermitian);
    Index<SymmGroup> right_i = ket_tensor.col_dim(),
                     out_left_i = physical_i * left_i;
    Index<SymmGroup> right_i_bra = bra_tensor.col_dim();
    common_subset(out_left_i, right_i);
    ProductBasis<SymmGroup> out_left_pb(physical_i, left_i);
    ProductBasis<SymmGroup> in_right_pb(physical_i, right_i,
                            boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                    -boost::lambda::_1, boost::lambda::_2));
    // Prepares output tensor
    MPSTensor<Matrix, SymmGroup> ret;
    ret.phys_i = bra_tensor.site_dim();
    ret.left_i = bra_tensor.row_dim();
    ret.right_i = bra_tensor.col_dim();
    auto loop_max = mpo.col_dim();
    bra_tensor.make_right_paired();
#ifdef USE_AMBIENT
        {
            block_matrix<Matrix, SymmGroup> empty;
            swap(ket_tensor.data(), empty); // deallocating mpstensor before exiting the stack
        }
        parallel::sync();
        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, left.aux_dim(), mpo.col_dim());

        parallel_for(index_type b2, parallel::range<index_type>(0,loop_max), {
            abelian::lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.data().basis(), bra_tensor.data().basis(), right_i, out_left_i, in_right_pb, out_left_pb);
        });
        omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
            contr_grid.multiply_column(b2, right[b2]);
        });
        t.clear();
        parallel::sync();

        swap(ret.data(), contr_grid.reduce());
        parallel::sync();
#else
    omp_for(index_type b2, parallel::range<std::size_t>(0,loop_max), {
        ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);
        abelian::lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.data().basis(), bra_tensor.data().basis(), right_i, out_left_i, in_right_pb, out_left_pb,
                             isHermitian);
        block_matrix<Matrix, SymmGroup> tmp;
        if (mpo.herm_info.right_skip(b2) && isHermitian)
            gemm(contr_grid(0,0), adjoint(right[mpo.herm_info.right_conj(b2)]), tmp);
        else
            gemm(contr_grid(0,0), right[b2], tmp);
        contr_grid(0,0).clear();
        parallel_critical
        for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
            ret.data().match_and_add_block(tmp[k], tmp.basis().left_charge(k), tmp.basis().right_charge(k));
    });
#endif
    return ret;
}

} // namespace contraction

#endif
