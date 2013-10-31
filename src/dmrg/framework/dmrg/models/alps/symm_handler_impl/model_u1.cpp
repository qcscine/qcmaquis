/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/alps/symm_handler.hpp"
#include "dmrg/block_matrix/symmetry/u1.h"
	
// Symmetry dependent implementation

// U1 Symmetry
template <>
U1::charge init_charge<U1> (const alps::Parameters& parms, std::map<std::string, int> const& all_conserved_qn)
{
    typedef std::map<std::string, int> qn_map_type;
    assert(all_conserved_qn.size() == 1);
    qn_map_type::const_iterator it = all_conserved_qn.begin();
    alps::half_integer<short> tmp = alps::evaluate<double>(static_cast<std::string>(parms[it->first+"_total"]), parms);
    return U1::charge( detail::to_integer(tmp) );
}

template <>
U1::charge state_to_charge<U1>(alps::site_state<short> const & state, alps::SiteBasisDescriptor<short> const& b,
                               std::map<std::string, int> const& all_conserved_qn)
{
    typedef std::map<std::string, int> qn_map_type;
    U1::charge c = U1::IdentityCharge;
    for (typename alps::SiteBasisDescriptor<short>::const_iterator it = b.begin(); it != b.end(); ++it) {
        qn_map_type::const_iterator match = all_conserved_qn.find(it->name());
        if (match != all_conserved_qn.end())
            c = detail::to_integer( get_quantumnumber(state, it->name(), b) );
    }
    return c;
}
