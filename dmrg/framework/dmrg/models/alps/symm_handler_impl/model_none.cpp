/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/alps/symm_handler.hpp"
#include "dmrg/block_matrix/symmetry/none.h"

// Symmetry dependent implementation

// TrivialGroup Symmetry
template <>
TrivialGroup::charge init_charge<TrivialGroup> (const alps::Parameters& parms, std::map<std::string, int> const& all_conserved_qn)
{
    return TrivialGroup::IdentityCharge;
}

template <>
TrivialGroup::charge state_to_charge<TrivialGroup>(alps::site_state<short> const & state, alps::SiteBasisDescriptor<short> const& b,
                                                   std::map<std::string, int> const& all_conserved_qn)
{
    return TrivialGroup::IdentityCharge;
}
