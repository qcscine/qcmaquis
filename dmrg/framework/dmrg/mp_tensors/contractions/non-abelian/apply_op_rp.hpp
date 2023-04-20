/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_SU2_APPLY_OP_RP_HPP
#define CONTRACTIONS_SU2_APPLY_OP_RP_HPP

#include "dmrg/block_matrix/symmetry/gsl_coupling.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/contractions/non-abelian/functors.h"

namespace contraction {
namespace SU2 {

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel_rp(size_t b2,
                        ContractionGrid<Matrix, SymmGroup>& contr_grid,
                        Boundary<OtherMatrix, SymmGroup> const & left,
                        BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> const & left_mult_mps,
                        MPOTensor<Matrix, SymmGroup> const & mpo,
                        DualIndex<SymmGroup> const & ket_basis,
                        Index<SymmGroup> const & right_i,
                        Index<SymmGroup> const & out_left_i,
                        ProductBasis<SymmGroup> const & in_right_pb,
                        ProductBasis<SymmGroup> const & out_left_pb)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename SymmGroup::charge charge;

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();

            block_matrix<Matrix, SymmGroup> local;
            block_matrix<Matrix, SymmGroup> const & T = left_mult_mps.at(b1, local);

            MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo.at(b1,b2);

        for (std::size_t op_index = 0; op_index < access.size(); ++op_index)
        {
            typename operator_selector<Matrix, SymmGroup>::type const & W = access.op(op_index);
            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);

            int a = mpo.left_spin(b1).get(), k = W.spin().get(), ap = mpo.right_spin(b2).get();

            for (size_t t_block = 0; t_block < T.n_blocks(); ++t_block) {

                charge lc = T.basis().left_charge(t_block);
                charge rc = T.basis().right_charge(t_block);

                const_iterator it = ket_basis.left_lower_bound(rc); //ket_basis comes transposed!
                charge mc = it->rc;

                for (size_t w_block = 0; w_block < W.basis().size(); ++w_block)
                {
                    charge phys_in = W.basis().left_charge(w_block);
                    charge phys_out = W.basis().right_charge(w_block);

                    charge out_r_charge = SymmGroup::fuse(rc, phys_in);
                    size_t rb = right_i.position(out_r_charge);
                    if (rb == right_i.size()) continue;

                    charge out_l_charge = SymmGroup::fuse(lc, phys_out);
                    if (!::SU2::triangle(SymmGroup::spin(out_r_charge), ap, SymmGroup::spin(out_l_charge))) continue;
                    if (!right_i.has(out_l_charge)) continue; // can also probe out_left_i, but right_i has the same charges

                    charge out_r_charge_rp = SymmGroup::fuse(out_r_charge, -phys_out);

                    size_t r_size = right_i[rb].second;

                    size_t o = ret.find_block(lc, out_r_charge_rp);
                    if ( o == ret.n_blocks() ) {
                        o = ret.insert_block(Matrix(T.basis().left_size(t_block),
                                             in_right_pb.size(-phys_out, out_r_charge)), lc, out_r_charge_rp);
                    }

                    int i = SymmGroup::spin(lc), ip = SymmGroup::spin(out_l_charge);
                    int j = SymmGroup::spin(mc), jp = SymmGroup::spin(out_r_charge);
                    int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);

                    typename Matrix::value_type couplings[4];
                    ::SU2::set_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip, access.scale(op_index), couplings);

                    size_t in_right_offset = in_right_pb(phys_in, out_r_charge);
                    size_t out_right_offset = in_right_pb(phys_out, out_r_charge);
                    size_t l_size = T.basis().left_size(t_block);

                    detail::rbtm_blocked<Matrix, SymmGroup>(T[t_block], ret[o], W, in_right_offset, out_right_offset, l_size, r_size, w_block, couplings);
                } // wblock
            } // lblock
        } // op_index
        } // b1
    }

} // namespace SU2
} // namespace contraction

#endif
