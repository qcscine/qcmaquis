/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPS_SECTORS_H
#define MPS_SECTORS_H

#include "dmrg/block_matrix/indexing.h"

template<class T>
T tri_min(T a, T b, T c)
{
    return std::min(std::min(a, b),
                    std::min(a, c));
}

template <class SymmGroup>
std::vector<Index<SymmGroup> > allowed_sectors(std::size_t L,
                                               Index<SymmGroup> const& phys,
                                               typename SymmGroup::charge right_end,
                                               std::size_t Mmax)
{
    bool finitegroup = SymmGroup::finite;
    
    Index<SymmGroup> physc = phys;
    physc.sort();
    typename SymmGroup::charge cmax = physc.begin()->first;
    typename SymmGroup::charge cmin = physc.rbegin()->first;
    if (cmin > cmax) std::swap(cmin, cmax);
    
    typename SymmGroup::charge cmaxL=SymmGroup::IdentityCharge, cminL=SymmGroup::IdentityCharge;
    for (int i = 1; i < L; ++i) {
        cmaxL = SymmGroup::fuse(cmaxL, cmax);
        cminL = SymmGroup::fuse(cminL, cmin);
    }
    
    Index<SymmGroup> l_triv, r_triv;
    l_triv.insert( std::make_pair(SymmGroup::IdentityCharge, 1) );
    r_triv.insert( std::make_pair(right_end, 1) );
    
    std::vector<Index<SymmGroup> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
    left_allowed[0] = l_triv;
    right_allowed[L] = r_triv;
    
    typename SymmGroup::charge cmaxi=cmaxL, cmini=cminL;
    for (int i = 1; i < L+1; ++i) {
        left_allowed[i] = phys * left_allowed[i-1];
        typename Index<SymmGroup>::iterator it = left_allowed[i].begin();
        while ( it != left_allowed[i].end() )
        {
            if (!finitegroup && SymmGroup::fuse(it->first, cmaxi) < right_end)
                it = left_allowed[i].erase(it);
            else if (!finitegroup && SymmGroup::fuse(it->first, cmini) > right_end)
                it = left_allowed[i].erase(it);
            else {
                it->second = std::min(Mmax, it->second);
                ++it;
            }
        }
        cmaxi = SymmGroup::fuse(cmaxi, -cmax);
        cmini = SymmGroup::fuse(cmini, -cmin);
    }
    cmaxi=cmaxL; cmini=cminL;
    for (int i = L-1; i >= 0; --i) {
        right_allowed[i] = adjoin(phys) * right_allowed[i+1];
        
        typename Index<SymmGroup>::iterator it = right_allowed[i].begin();
        while ( it != right_allowed[i].end() )
        {
            if (!finitegroup && SymmGroup::fuse(it->first, -cmaxi) > SymmGroup::IdentityCharge)
                it = right_allowed[i].erase(it);
            else if (!finitegroup && SymmGroup::fuse(it->first, -cmini) < SymmGroup::IdentityCharge)
                it = right_allowed[i].erase(it);
            else {
                it->second = std::min(Mmax, it->second);
                ++it;
            }
        }
        cmaxi = SymmGroup::fuse(cmaxi, -cmax);
        cmini = SymmGroup::fuse(cmini, -cmin);
        
    }
    
    for (int i = 0; i < L+1; ++i) {
        allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
        for (typename Index<SymmGroup>::iterator it = allowed[i].begin();
             it != allowed[i].end(); ++it)
            it->second = tri_min(Mmax,
                                 left_allowed[i].size_of_block(it->first),
                                 right_allowed[i].size_of_block(it->first));
    }
    
    return allowed;
}


#endif
