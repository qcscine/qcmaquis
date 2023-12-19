/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/alps/symm_handler.hpp"
#include "dmrg/block_matrix/symmetry/nu1pg.h"

// Symmetry dependent implementation

// TwoU1PG Symmetry
template <>
TwoU1PG::charge state_to_charge<TwoU1PG>(alps::site_state<short> const & state, alps::SiteBasisDescriptor<short> const& b,
                                     std::map<std::string, int> const& all_conserved_qn)
{
    typedef std::map<std::string, int> qn_map_type;
    TwoU1PG::charge c = TwoU1PG::IdentityCharge;
    for (alps::SiteBasisDescriptor<short>::const_iterator it = b.begin(); it != b.end(); ++it) {
        qn_map_type::const_iterator match = all_conserved_qn.find(it->name());
        if (match != all_conserved_qn.end())
            c[match->second] = detail::to_integer( get_quantumnumber(state, it->name(), b) );
    }
    return c;
}

template <>
TwoU1PG::charge init_charge<TwoU1PG> (const alps::Parameters& parms, std::map<std::string, int> const& all_conserved_qn)
{
    typedef std::map<std::string, int> qn_map_type;
    assert(all_conserved_qn.size() == 2);

    TwoU1PG::charge c = TwoU1PG::IdentityCharge;
    for (qn_map_type::const_iterator it=all_conserved_qn.begin(); it!=all_conserved_qn.end(); ++it) {
        alps::half_integer<short> tmp = alps::evaluate<double>(static_cast<std::string>(parms[it->first+"_total"]), parms);
        c[it->second] = detail::to_integer(tmp);
    }
    
    return c;
}
