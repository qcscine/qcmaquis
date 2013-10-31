/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/alps/model_symm.hpp"
#include "dmrg/block_matrix/symmetry/z2.h"
	
// Symmetry dependent implementation

// Z2 Symmetry
template <>
Ztwo::charge init_charge<Ztwo> (const alps::Parameters& parms, std::map<std::string, int> const& all_conserved_qn)
{
    assert(all_conserved_qn.size() == 1);
    qn_map_type::const_iterator it = all_conserved_qn.begin();
    int tmp = alps::evaluate<double>(static_cast<std::string>(parms[it->first+"_total"]), parms);
    if (!(tmp == 0 || tmp == 1))
        throw std::runtime_error("Invalid value for " + it->first + "_total");
    return (tmp == 1) ? Ztwo::Minus : Ztwo::Plus;
}

template <>
Ztwo::charge state_to_charge<Ztwo>(alps::site_state<short> const & state, alps::SiteBasisDescriptor<short> const& b,
                               std::map<std::string, int> const& all_conserved_qn)
{
    typedef std::map<std::string, int> qn_map_type;
    int tmp = 0;
    for (typename alps::SiteBasisDescriptor<short>::const_iterator it = b.begin(); it != b.end(); ++it) {
        qn_map_type::const_iterator match = all_conserved_qn.find(it->name());
        if (match != all_conserved_qn.end())
            tmp = detail::to_integer( get_quantumnumber(state, it->name(), b) );
    }
    return (tmp % 2 == 0 ? Ztwo::Plus : Ztwo::Minus);
}
