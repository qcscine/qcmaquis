/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_STATE_MPS_H
#define MAQUIS_DMRG_STATE_MPS_H

#include "dmrg/mp_tensors/mps.h"

#include <boost/tuple/tuple.hpp>

template <class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> state_mps(std::vector<boost::tuple<typename SymmGroup::charge, size_t> > const & state,
                                 Index<SymmGroup> const& phys)
{
    typedef typename SymmGroup::charge charge;
    typedef boost::tuple<charge, size_t> local_state;
    
    MPS<Matrix, SymmGroup> mps(state.size());
    
    Index<SymmGroup> curr_i;
    curr_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
    size_t curr_b = 0;
    for (int i=0; i<state.size(); ++i)
    {
        charge newc = SymmGroup::fuse(curr_i[0].first, boost::get<0>(state[i]));
        size_t news = 1;
        Index<SymmGroup> new_i;
        new_i.insert(std::make_pair(newc, news));
        ProductBasis<SymmGroup> left(phys, curr_i);
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys, curr_i, new_i, false, 0);
        size_t b_in = left(boost::get<0>(state[i]), curr_i[0].first) + boost::get<1>(state[i]) * curr_i[0].second + curr_b;
        size_t b_out = 0;
        
        mps[i].make_left_paired();
        block_matrix<Matrix, SymmGroup> & block = mps[i].data();
        Matrix & m = block(SymmGroup::fuse(curr_i[0].first, boost::get<0>(state[i])), new_i[0].first);
        m(b_in, b_out) = 1.;
        
        curr_i = new_i;
        curr_b = b_out;
    }
    return mps;
}


#endif
