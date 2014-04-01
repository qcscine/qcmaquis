/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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
