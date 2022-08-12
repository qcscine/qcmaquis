/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_STATE_MPS_H
#define MAQUIS_DMRG_STATE_MPS_H

#include "dmrg/mp_tensors/mps.h"
#include "mps_sectors.h"
#include <boost/tuple/tuple.hpp>

template <class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> state_mps(std::vector<boost::tuple<typename SymmGroup::charge, size_t> > const & state,
                                 std::vector<Index<SymmGroup> > const& phys_dims,
                                 std::vector<int> const& site_type,
                                 typename SymmGroup::charge right_end = SymmGroup::IdentityCharge,
                                 int mdim=1)
{
    // Types and variable definition
    typedef typename SymmGroup::charge charge;
    MPS<Matrix, SymmGroup> mps(state.size());
    Index<SymmGroup> curr_i;
    std::vector< Index<SymmGroup> > allowed = allowed_sectors(site_type, phys_dims, right_end, mdim);
    // -- MAIN LOOP --
    // The overall structure of the algorithm is the following:
    // One first generates the index (IdentityCharge,mdim), where mdim is the number of renormalized
    // block states (1 in the default case). This is the index with which we start from the left.
    // For a given left index and a given physical basis state, the "acceptable" QN for the right 
    // renormalized basis are univocally determined by the product of the two set of QNs.
    // Therefore, since we have a single ONV as a starting point, we will have
    // all 1x1 tensors, but the physical dimensions will be in general != from 0.
    curr_i.insert(std::make_pair(SymmGroup::IdentityCharge, mdim));
    for (int i = 0; i < state.size(); ++i)
    {
        // Computes the symmetry block of the next dimension and HARDCODED its value to 1
        charge newc = SymmGroup::fuse(curr_i[0].first, boost::get<0>(state[i]));
        Index<SymmGroup> new_i;
        new_i.insert(std::make_pair(newc, mdim));
        // Get the product basis between the physical basis and the symmetry block of the left renormalized basis
        ProductBasis<SymmGroup> left(phys_dims[site_type[i]], allowed[i]);
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], false, 0);
        // Finds out where to put the 1.0 in the MPS. Retrieve, from the ProductBasis object, how the row index was
        // decomposed in terms of left auxiliary basis and physical basis.
        size_t b_in = left(boost::get<0>(state[i]), curr_i[0].first) + boost::get<1>(state[i]) * curr_i[0].second;
        assert (allowed[i+1].has(newc));
        size_t b_out = 0;
        mps[i].make_left_paired();
        // Populates the MPS
        block_matrix<Matrix, SymmGroup> & block = mps[i].data();
        Matrix &m = block(newc, new_i[0].first);
        m(b_in, b_out) = 1.;
        curr_i = new_i;
    }
    return mps;
}


#endif
