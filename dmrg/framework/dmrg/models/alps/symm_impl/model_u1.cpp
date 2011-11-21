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
    
    // U1 Symmetry
    template <>
    U1::charge init_charge<U1> (const alps::Parameters& parms, std::vector<std::pair<std::size_t, std::string> > const& qn)
    {
        assert(qn.size() == 1);
        U1::charge c = U1::IdentityCharge;
        if (parms.defined(qn[0].second+"_total")) {
            alps::half_integer<short> tmp = alps::evaluate<double>(static_cast<std::string>(parms[qn[0].second+"_total"]), parms);
            c = details::to_integer(tmp);
        }
        return c;
    }

    template <>
    U1::charge convert_alps<U1> (alps::site_state<short> const & state, std::vector<std::pair<std::size_t, std::string> > const& qn)
    {
        assert(qn.size() == 1);
        return details::to_integer( get_quantumnumber(state, qn[0].first) );
    }
    
    template <>
    std::map<U1::charge,std::size_t> init_qn_charges<U1>
    (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
    {
        assert(conserved_qn.size() == 1);
        std::map<U1::charge,std::size_t> ret;
        for (int i=0; i<states.size(); ++i)
            ret[convert_alps<U1>(states[i], conserved_qn)] = 1;
        return ret;
    }

    template <>
    std::map<alps::site_state<short>, std::pair<U1::charge, std::size_t> > init_coords<U1>
    (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
    {
        assert(conserved_qn.size() == 1);
        std::map<alps::site_state<short>, std::pair<U1::charge, std::size_t> > ret;
        for (std::size_t i=0; i<states.size(); ++i)
            ret[states[i]] = std::make_pair(convert_alps<U1>(states[i], conserved_qn), 0);
        return ret;
    }

}
