/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_IDENTITY_MPS_H
#define MAQUIS_DMRG_IDENTITY_MPS_H

#include "dmrg/block_matrix/grouped_symmetry.h"

template <class Matrix, class InSymm>
MPS<Matrix, typename grouped_symmetry<InSymm>::type> identity_dm_mps(std::size_t L, Index<InSymm> const& phys_psi,
                                                                     std::vector<Index<typename grouped_symmetry<InSymm>::type> > const& allowed)
{
    MPOTensor<Matrix, InSymm> t(1,1);
    t.set(0, 0, identity_matrix<typename operator_selector<Matrix, SymmGroup>::type>(phys_psi));

    MPO<Matrix, InSymm> mpo(L, t);
    return mpo_to_smps_group(mpo, phys_psi, allowed);
}

#endif
