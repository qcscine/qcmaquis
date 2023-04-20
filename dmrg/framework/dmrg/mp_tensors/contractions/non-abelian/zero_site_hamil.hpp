/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_SU2_ZERO_SITE_HAMIL_HPP
#define CONTRACTIONS_SU2_ZERO_SITE_HAMIL_HPP

namespace contraction {

template<class Matrix, class OtherMatrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
zerosite_hamil_lbtm(block_matrix<Matrix, SymmGroup> bra_tensor, block_matrix<Matrix, SymmGroup> ket_tensor, Boundary<OtherMatrix, SymmGroup> const & left, 
                    Boundary<OtherMatrix, SymmGroup> const & right, MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right,
                    bool isHermitian);

template<class Matrix, class OtherMatrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
zerosite_hamil_lbtm(block_matrix<Matrix, SymmGroup> ket_tensor, Boundary<OtherMatrix, SymmGroup> const & left, 
                    Boundary<OtherMatrix, SymmGroup> const & right, MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right,
                    bool isHermitian);

template<class Matrix, class OtherMatrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
Engine<Matrix, OtherMatrix, SymmGroup, typename symm_traits::enable_if_su2_t<SymmGroup>>::
zerosite_hamil2(block_matrix<Matrix, SymmGroup> bra_tensor, block_matrix<Matrix, SymmGroup> ket_tensor, Boundary<OtherMatrix, SymmGroup> const & left,
                Boundary<OtherMatrix, SymmGroup> const & right, MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right,
                bool isHermitian)
{
    return zerosite_hamil_lbtm(bra_tensor, ket_tensor, left, right, mpo_left, mpo_right, isHermitian);
}

template<class Matrix, class OtherMatrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
Engine<Matrix, OtherMatrix, SymmGroup, typename symm_traits::enable_if_su2_t<SymmGroup>>::
zerosite_hamil2(block_matrix<Matrix, SymmGroup> ket_tensor, Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right, bool isHermitian)
{
    return zerosite_hamil_lbtm(ket_tensor, left, right, mpo_left, mpo_right, isHermitian);
}

template<class Matrix, class OtherMatrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
zerosite_hamil_lbtm_kernel(block_matrix<Matrix, SymmGroup> const & bra_tensor, block_matrix<Matrix, SymmGroup> ket_tensor,
                           Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                           MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right,
                           bool isHermitian)
{
    // General types and variables
    typedef typename SymmGroup::charge charge;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    typedef typename Matrix::value_type value_type;
    block_matrix<Matrix, SymmGroup> ret;
    // Contraction with the boundary
    contraction::common::BoundaryMPSProduct<Matrix, OtherMatrix, SymmGroup, ::SU2::SU2Gemms> t(ket_tensor, left, mpo_right, isHermitian, false);
    index_type loop_max = mpo_right.row_dim();
    // Final contraction with the right boundary
    omp_for(index_type b2, parallel::range<index_type>(0,loop_max), {
        // Variables and types
        block_matrix<Matrix, SymmGroup> tmp, tmp2, local;
        // Contraction
        if (mpo_right.herm_info.left_skip(b2) && isHermitian) {
            //std::vector<value_type> phases = ::contraction::common::conjugate_phases(adjoint(right[mpo_left.herm_info.right_conj(b2)]), mpo_left, b2, false, true);
            //::SU2::gemm_trim(t.at(b2, local), adjoint(right[mpo_left.herm_info.right_conj(b2)]), tmp, phases, false);
            ::SU2::gemm_trim(t.at(b2, local), adjoint(right[mpo_left.herm_info.right_conj(b2)]), tmp, std::vector<value_type>(adjoint(right[mpo_left.herm_info.right_conj(b2)]).n_blocks(), 1.), false);
        }
        else {
            ::SU2::gemm_trim(t.at(b2, local), right[b2], tmp, std::vector<value_type>(t.at(b2, local).n_blocks(), 1.), true);
        }
        parallel_critical
        for (std::size_t k = 0; k < tmp.n_blocks(); ++k)
            ret.match_and_add_block(tmp[k], tmp.basis().left_charge(k), tmp.basis().right_charge(k));
    });
    return ret;
}

template<class Matrix, class OtherMatrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
zerosite_hamil_lbtm(block_matrix<Matrix, SymmGroup> ket_tensor, Boundary<OtherMatrix, SymmGroup> const & left,
                     Boundary<OtherMatrix, SymmGroup> const & right, MPOTensor<Matrix, SymmGroup> const & mpo_left,
                     MPOTensor<Matrix, SymmGroup> const & mpo_right, bool isHermitian)
{
    return zerosite_hamil_lbtm_kernel<Matrix, OtherMatrix, SymmGroup>(ket_tensor, ket_tensor, left, right, mpo_left, mpo_right, isHermitian);
}

template<class Matrix, class OtherMatrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
zerosite_hamil_lbtm(block_matrix<Matrix, SymmGroup> bra_tensor, block_matrix<Matrix, SymmGroup> ket_tensor,
                    Boundary<OtherMatrix, SymmGroup> const & left, Boundary<OtherMatrix, SymmGroup> const & right,
                    MPOTensor<Matrix, SymmGroup> const & mpo_left, MPOTensor<Matrix, SymmGroup> const & mpo_right,
                    bool isHermitian)
{
    //return zerosite_hamil_lbtm_kernel<Matrix, OtherMatrix, SymmGroup>(bra_tensor, ket_tensor, left, right, mpo_left, mpo_right, isHermitian);
    throw std::runtime_error("Zero site Hamiltonian for bra != ket NYI");
}

} // namespace contraction

#endif
