/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef JOINS_H
#define JOINS_H

#include "dmrg/mp_tensors/mpstensor.h"

namespace MPSJoin {

// Forward declaraation of ADdMPSAsBlock
template<class Matrix, class SymmGroup>
void AddMPSAsBlock(MPSTensor<Matrix, SymmGroup> const& input, MPSTensor<Matrix, SymmGroup> const& ref_for_offset,
                   MPSTensor<Matrix, SymmGroup>& output, ProductBasis<SymmGroup> const& out_left_pb,
                   bool do_left_offset=false, bool do_right_offset=false);

/*!
 * Routine to sum two MPSs. Note that the sum of two MPSs is given by the direct sum of the MPSs (i.e., a block-shaped
 * sum) and not simply by the sum of the single components.
 *
 * \tparam Matrix: type of the two-dim matrix underlying the MPS implementation.
 * \tparam SymmGroup: class representing the symmetry group of the molecule.
 * \param m1: first matrix to be summed.
 * \param m2: second matrix to be summed.
 * \param boundary_f: enum used to check if the MPSTensor to join is at the beginning or at the end of the chain, in
 *                    which case no offset must be applied.
 * \return MPS representing m1 + m2.
 */
template <class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> join(MPSTensor<Matrix, SymmGroup> const & m1, MPSTensor<Matrix, SymmGroup> const & m2,
                                  boundary_flag_t boundary_f=no_boundary_f)
{
    // If in one of the two MPSs one of the physics index is missing, add it.
    MPSTensor<Matrix, SymmGroup> ret;
    auto phys_i = m1.site_dim();
    if (m1.site_dim() != m2.site_dim())
        for (auto it = m2.site_dim().begin(); it != m2.site_dim().end(); ++it)
            if (!phys_i.has(it->first))
               phys_i.insert(*it);  
    // MPS are made left paired before merging them.
    m1.make_left_paired();
    m2.make_left_paired();
    // -- INDEXES MANAGEMENT --
    // The sum changes only the auxiliary indexes, not the physical ones that are not altered.
    ret.phys_i = phys_i;
    // Computes the new value of the dimension for the left renormalized basis as the sum of the dimensions of the
    // basis in m1 and m2.
    ret.left_i = m1.left_i;
    if (boundary_f != l_boundary_f) {
        for (auto it = m2.left_i.begin(); it != m2.left_i.end(); ++it)
            if (ret.left_i.has(it->first))
                ret.left_i[ret.left_i.position(it->first)].second += it->second;
            else
                ret.left_i.insert(*it);
    }
    // Does the same for the right index.
    ret.right_i = m1.right_i;
    if (boundary_f != r_boundary_f) {
        for (auto it = m2.right_i.begin(); it != m2.right_i.end(); ++it)
            if (ret.right_i.has(it->first))
                ret.right_i[ret.right_i.position(it->first)].second += it->second;
            else
                ret.right_i.insert(*it);
    }
    // Generates the dimension of the summed MPS tensor.
    ProductBasis<SymmGroup> out_left_pb(phys_i, ret.left_i);
    Index<SymmGroup> const& out_right = ret.right_i;
    using std::size_t;
    // -- SUM OF THE BLOCKS OF THE MPSTENSORS --
    MPSJoin::AddMPSAsBlock(m1, m1, ret, out_left_pb);
    MPSJoin::AddMPSAsBlock(m2, m1, ret, out_left_pb, boundary_f != l_boundary_f, boundary_f != r_boundary_f);
    // check right_pairing
    assert( weak_equal(ret.right_i, ret.data().right_basis()) );
    return ret;
}

/*!
 * Routine to merge three MPSs (M1, M12, M2) in a single super-MPS according to the following structure:
 * +-         -+
 * |           |
 * | M1     0  |
 * |           |     ==> This is required in the calculation of the tangent-space vector.
 * | M12   M2  |
 * |           |
 * +--        -+
 *
 * \tparam Matrix: type of the two-dim matrix underlying the MPS implementation.
 * \tparam SymmGroup: class representing the symmetry group of the molecule.
 * \param m1: first matrix to be summed.
 * \param m2: second matrix to be summed.
 * \param m12: off-diagonal block.
 * \return MPS representing the super-MPS.
 *
 * Note that the enum for the initial/final site is not anymore needed since it is not possible to have more than
 * two blocks at the extremal points of the DMRG lattice.
 */
template <class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> join(MPSTensor<Matrix, SymmGroup> const & m1, MPSTensor<Matrix, SymmGroup> const & m2,
                                  MPSTensor<Matrix, SymmGroup> const & m12)
{
    // If in one of the two MPSs one of the physics index is missing, add it.
    MPSTensor<Matrix, SymmGroup> ret;
    auto phys_i = m1.site_dim();
    if (m1.site_dim() != m2.site_dim())
        for (auto it = m2.site_dim().begin(); it != m2.site_dim().end(); ++it)
            if (!phys_i.has(it->first))
                phys_i.insert(*it);
    //
    if (m1.site_dim() != m12.site_dim())
        for (auto it = m12.site_dim().begin(); it != m12.site_dim().end(); ++it)
            if (!phys_i.has(it->first))
                phys_i.insert(*it);
    // MPS are made left paired before merging them.
    m1.make_left_paired();
    m2.make_left_paired();
    m12.make_left_paired();
    // Now the coherence of the dimensions of the m12 tensor with the other two tensors
    for (auto it = m12.left_i.begin(); it != m12.left_i.end(); it++)
        if (m2.left_i.has(it->second)) {
            if (m2.left_i.size_of_block(it->second) != m12.left_i.size_of_block(it->second))
                throw std::runtime_error("Dimension of symmetry blocks not coherent.");
        } else {
            throw std::runtime_error("Symmetry block not found.");
        }
    //
    for (auto it = m12.left_i.begin(); it != m12.left_i.end(); it++)
        if (m1.right_i.has(it->second)) {
            if (m2.right_i.size_of_block(it->second) != m12.right_i.size_of_block(it->second))
                throw std::runtime_error("Dimension of symmetry blocks not coherent.");
        } else {
            throw std::runtime_error("Symmetry block not found.");
        }
    // -- INDEXES MANAGEMENT --
    // The sum changes only the auxiliary indexes, not the physical ones that are not altered.
    ret.phys_i = phys_i;
    // Computes the new value of the dimension for the left renormalized basis as the sum of the dimensions of the
    // basis in m1 and m2.
    ret.left_i = m1.left_i;
    for (auto it = m2.left_i.begin(); it != m2.left_i.end(); ++it)
        if (ret.left_i.has(it->first))
            ret.left_i[ret.left_i.position(it->first)].second += it->second;
        else
            ret.left_i.insert(*it);
    // Does the same for the right index.
    ret.right_i = m1.right_i;
    for (auto it = m2.right_i.begin(); it != m2.right_i.end(); ++it)
        if (ret.right_i.has(it->first))
            ret.right_i[ret.right_i.position(it->first)].second += it->second;
        else
            ret.right_i.insert(*it);
    // Generates the dimension of the summed MPS tensor.
    ProductBasis<SymmGroup> out_left_pb(phys_i, ret.left_i);
    Index<SymmGroup> const& out_right = ret.right_i;
    using std::size_t;
    // Add the various blocks.
    MPSJoin::AddMPSAsBlock(m1, m1, ret, out_left_pb);
    MPSJoin::AddMPSAsBlock(m2, m1, ret, out_left_pb, true, true);
    MPSJoin::AddMPSAsBlock(m12, m1, ret, out_left_pb, true, false);
    // check right_pairing
    assert( weak_equal(ret.right_i, ret.data().right_basis()) );
    return ret;
}

/*!
 * Wrapper to add a block to a given super-block adding an offset based on a reference MPS.
 * Used to avoid code repetition in the part of the code above.
 * \tparam Matrix: type of the two-dim matrix underlying the MPS implementation.
 * \tparam SymmGroup: class representing the symmetry group of the molecule.
 * \param input: MPSTensor to be added.
 * \param ref_for_offset: MPSTensor taken as reference to calculate the offsets (see M1 in [join])
 * \param output: output MPSTensor
 * \param out_left_pb: product basis for the output MPS.
 * \param do_left_offset: apply left offset.
 * \param do_right_offset: apply right offset.
 */

template<class Matrix, class SymmGroup>
void AddMPSAsBlock(MPSTensor<Matrix, SymmGroup> const& input, MPSTensor<Matrix, SymmGroup> const& ref_for_offset,
                   MPSTensor<Matrix, SymmGroup>& output, ProductBasis<SymmGroup> const& out_left_pb,
                   bool do_left_offset, bool do_right_offset)
{
    // Input parameters
    Index<SymmGroup> const& phys_i_m = input.site_dim();
    ProductBasis<SymmGroup> in_left(phys_i_m, input.row_dim());
    // Loop over all symmetry blocks. Remember that both matrices are left paired.
    for (size_t b = 0; b < input.data().n_blocks(); ++b) {
        // Retrieves the position of the symmetry block in the output tensor.
        typename SymmGroup::charge const& sl_charge = input.data().basis().left_charge(b); // phys + left
        typename SymmGroup::charge const& r_charge = input.data().basis().right_charge(b); // right
        size_t out_r_charge_i = output.right_i.position(r_charge);
        if (!output.data().has_block(sl_charge, r_charge))
            output.data().insert_block(Matrix(out_left_pb.size(sl_charge),
                                              output.right_i[out_r_charge_i].second),
                                              sl_charge, r_charge);
        Matrix& nb = output.data()(sl_charge, r_charge);
        size_t in_r_size = input.data().basis().right_size(b);
        size_t out_r_offset = 0;
        // Calculates the offset, starting from which the matrix is populated
        if (do_right_offset)
            out_r_offset += ref_for_offset.col_dim().size_of_block(r_charge, true);
        // Now copies the actual data. Remember that a given block, in the left pairing case, is the sum of the
        // contribution of different (phys_i, left_i) pairs.
        for (size_t s = 0; s < phys_i_m.size(); ++s) {
            typename SymmGroup::charge const& s_charge = phys_i_m[s].first;
            typename SymmGroup::charge l_charge = SymmGroup::fuse(sl_charge, -s_charge); // left
            if (!input.row_dim().has(l_charge))
                continue;
            auto in_l_size = input.row_dim().size_of_block(l_charge, true);
            auto in_l_offset = in_left(s_charge, l_charge);
            auto out_l_size = output.row_dim().size_of_block(l_charge, true);
            auto out_l_offset = out_left_pb(s_charge, l_charge);
            // Second offset
            if (do_left_offset)
                out_l_offset += ref_for_offset.row_dim().size_of_block(l_charge, true);
            for (size_t ss = 0; ss < output.phys_i[s].second; ++ss)
                copy_block(input.data()[b], in_l_offset + ss * in_l_size, 0, nb, out_l_offset + ss * out_l_size,
                           out_r_offset, in_l_size, in_r_size);
        }
    }
}

}


#endif
