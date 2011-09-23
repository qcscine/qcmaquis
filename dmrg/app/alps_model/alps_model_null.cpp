/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "app/alps_model.hpp"

namespace app {
	
	// Symmetry dependent implementation
    
    // NullGroup Symmetry
    template <>
    NullGroup::charge init_charge<NullGroup> (const alps::Parameters& parms, std::vector<std::pair<std::size_t, std::string> > const& qn)
    {
        return NullGroup::SingletCharge;
    }
    
    template <>
    NullGroup::charge convert_alps<NullGroup> (alps::site_state<short> const & state, std::vector<std::pair<std::size_t, std::string> > const& qn)
    {
        return NullGroup::SingletCharge;
    }

    template <>
    std::map<NullGroup::charge,std::size_t> init_qn_charges<NullGroup>
    (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
    {
        std::map<NullGroup::charge, std::size_t> ret;
        ret[convert_alps<NullGroup>(states[0], conserved_qn)] = states.size();
        return ret;
    }
    
    template <>
    std::map<alps::site_state<short>, std::pair<NullGroup::charge, std::size_t> > init_coords<NullGroup>
    (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
    {
        std::map<alps::site_state<short>, std::pair<NullGroup::charge, std::size_t> > ret;
        for (std::size_t i=0; i<states.size(); ++i)
            ret[states[i]] = std::make_pair(NullGroup::SingletCharge, i);
        return ret;
    }
	
}
