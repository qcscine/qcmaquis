/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/alps/model.hpp"

namespace app {
	
	// Symmetry dependent implementation
    
    // TrivialGroup Symmetry
    template <>
    TrivialGroup::charge init_charge<TrivialGroup> (const alps::Parameters& parms, std::vector<std::pair<std::size_t, std::string> > const& qn)
    {
        return TrivialGroup::IdentityCharge;
    }
    
    template <>
    TrivialGroup::charge convert_alps<TrivialGroup> (alps::site_state<short> const & state, std::vector<std::pair<std::size_t, std::string> > const& qn)
    {
        return TrivialGroup::IdentityCharge;
    }

    template <>
    std::map<TrivialGroup::charge,std::size_t> init_qn_charges<TrivialGroup>
    (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
    {
        std::map<TrivialGroup::charge, std::size_t> ret;
        ret[convert_alps<TrivialGroup>(states[0], conserved_qn)] = states.size();
        return ret;
    }
    
    template <>
    std::map<alps::site_state<short>, std::pair<TrivialGroup::charge, std::size_t> > init_coords<TrivialGroup>
    (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
    {
        std::map<alps::site_state<short>, std::pair<TrivialGroup::charge, std::size_t> > ret;
        for (std::size_t i=0; i<states.size(); ++i)
            ret[states[i]] = std::make_pair(TrivialGroup::IdentityCharge, i);
        return ret;
    }
	
}
