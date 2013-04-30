/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_BLOCK_MATRIX_GROUPED_SYMMETRY_H
#define MAQUIS_DMRG_BLOCK_MATRIX_GROUPED_SYMMETRY_H

#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/indexing.h"

//// TRAITS

template <class SymmGroup>
struct grouped_symmetry;

template <>
struct grouped_symmetry<TrivialGroup> {
    typedef TrivialGroup type;
};

template <>
struct grouped_symmetry<U1> {
    typedef TwoU1 type;
};

//// GROUPPING FUNCTIONS

inline TrivialGroup::charge group(TrivialGroup::charge c1, TrivialGroup::charge c2)
{
    return TrivialGroup::IdentityCharge;
}

inline TwoU1::charge group(U1::charge c1, U1::charge c2)
{
    TwoU1::charge R;
    R[0] = c1; R[1] = c2;
    return R;
}

template<class SymmGroup>
Index<typename grouped_symmetry<SymmGroup>::type> group(Index<SymmGroup> const & i1,
                                             Index<SymmGroup> const & i2)
{
    typedef typename grouped_symmetry<SymmGroup>::type OutSymm;
    
    Index<OutSymm> ret;
    for (typename Index<SymmGroup>::const_iterator it1 = i1.begin(); it1 != i1.end(); ++it1)
        for (typename Index<SymmGroup>::const_iterator it2 = i2.begin(); it2 != i2.end(); ++it2)
        {
            ret.insert(std::make_pair(group(it1->first, it2->first), it1->second*it2->second));
        }
    return ret;
}


#endif
