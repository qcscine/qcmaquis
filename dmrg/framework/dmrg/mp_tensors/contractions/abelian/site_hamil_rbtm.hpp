/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef ABELIAN_SITE_HAMIL_RBTM
#define ABELIAN_SITE_HAMIL_RBTM

namespace contraction {

    template<class Matrix, class OtherMatrix, class SymmGroup, class SymmType>
    MPSTensor<Matrix, SymmGroup>
    Engine<Matrix, OtherMatrix, SymmGroup, SymmType>::
    site_hamil_rbtm(MPSTensor<Matrix, SymmGroup> ket_tensor,
                Boundary<OtherMatrix, SymmGroup> const & left,
                Boundary<OtherMatrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

        std::vector<block_matrix<Matrix, SymmGroup> > t
            = common::mps_times_boundary<Matrix, OtherMatrix, SymmGroup, Gemms>(ket_tensor, right, mpo);

        Index<SymmGroup> const & physical_i = ket_tensor.site_dim(),
                                 right_i = ket_tensor.col_dim();
        Index<SymmGroup> left_i = ket_tensor.row_dim(),
                         out_right_i = adjoin(physical_i) * right_i;

        common_subset(out_right_i, left_i);
        ProductBasis<SymmGroup> in_left_pb(physical_i, left_i);
        ProductBasis<SymmGroup> out_right_pb(physical_i, right_i,
                                             boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                 -boost::lambda::_1, boost::lambda::_2));
        block_matrix<Matrix, SymmGroup> collector;
        MPSTensor<Matrix, SymmGroup> ret;
        ret.phys_i = ket_tensor.site_dim(); ret.left_i = ket_tensor.row_dim(); ret.right_i = ket_tensor.col_dim();

        index_type loop_max = mpo.row_dim();
        omp_for(index_type b1, parallel::range<index_type>(0,loop_max), {

            block_matrix<Matrix, SymmGroup> tmp, tmp2;
            abelian::rbtm_kernel(b1, tmp, right, t, mpo, ket_tensor.data().basis(), left_i, out_right_i, in_left_pb, out_right_pb);

            gemm(transpose(left[b1]), tmp, tmp2);
            swap(tmp, tmp2);

            parallel_critical
            for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
                collector.match_and_add_block(tmp[k], tmp.basis().left_charge(k), tmp.basis().right_charge(k));
        });

        reshape_right_to_left_new(physical_i, left_i, right_i, collector, ret.data());
        return ret;
    }

} // namespace contraction

#endif
