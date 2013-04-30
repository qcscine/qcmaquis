/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_IDENTITY_MPS_H
#define MAQUIS_DMRG_IDENTITY_MPS_H

#include "dmrg/block_matrix/grouped_symmetry.h"

template <class Matrix, class InSymm>
MPS<Matrix, typename grouped_symmetry<InSymm>::type> identity_dm_mps(std::size_t L, Index<InSymm> const& phys_psi,
                                                                     std::vector<Index<typename grouped_symmetry<InSymm>::type> > const& allowed)
{
    MPOTensor<Matrix, InSymm> t(1,1);
    t(0,0) = identity_matrix<Matrix>(phys_psi);
    
    MPO<Matrix, InSymm> mpo(L, t);
    return mpo_to_smps_group(mpo, phys_psi, allowed);
}

#endif
