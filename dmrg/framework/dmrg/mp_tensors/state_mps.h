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

//
// STATE_MPS FUNCTION
// ------------------
//
// Function that takes in input:
// 1) a vector of tuples (Charge,Dimension), giving the state
// 2) a vector of indexes, giving the physical dimensions
// 3) a vector of integers, giving the site_type
// and returns in output the corresponding MPS
//

template <class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> state_mps(std::vector<boost::tuple<typename SymmGroup::charge, size_t> > const & state,
                                 std::vector<Index<SymmGroup> > const& phys_dims,
                                 std::vector<int> const& site_type,
                                 typename SymmGroup::charge right_end = SymmGroup::IdentityCharge,
                                 int const& mdim = 1)
{
    // Types definition
    typedef typename SymmGroup::charge charge;
    typedef boost::tuple<charge, size_t> local_state;
    MPS<Matrix, SymmGroup> mps(state.size());
    Index<SymmGroup> curr_i;
    std::vector< Index<SymmGroup> > allowed = allowed_sectors(site_type, phys_dims, right_end, mdim) ;
    // Generate the index (IdentityCharge,mdim), where mdim is the number of renormalized
    // block states (1 in the default case)
    curr_i.insert(std::make_pair(SymmGroup::IdentityCharge, mdim));
    size_t curr_b = 0;
    for (int i=0; i<state.size(); ++i)
    {
        charge newc = SymmGroup::fuse(curr_i[0].first, boost::get<0>(state[i]));
        size_t news = 1;
        Index<SymmGroup> new_i;
        new_i.insert(std::make_pair(newc, news));
        // Get the product basis between the physical basis and the symmetry block of the left renormalized basis
        ProductBasis<SymmGroup> left(phys_dims[site_type[i]], allowed[i]);
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], false, 0.);
        // Finds out where to put the 1.0 in the MPS. Retrieve, from the ProductBasis object, how the row index was
        // decomposed in terms of left, auxiliary basis and physical basis.
        size_t b_in = left(boost::get<0>(state[i]), curr_i[0].first) + boost::get<1>(state[i]) * curr_i[0].second + curr_b;
        size_t b_out = 0;
        
        mps[i].make_left_paired();
        block_matrix<Matrix, SymmGroup> & block = mps[i].data();
        Matrix &m = block(newc, new_i[0].first);
        m(b_in, b_out) = 1.;
        
        curr_i = new_i;
        curr_b = b_out;
    }
    return mps;
}


#endif
