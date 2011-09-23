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

    // TwoU1 Symmetry
    template <>
    TwoU1::charge convert_alps<TwoU1> (alps::site_state<short> const & state, std::vector<std::pair<std::size_t, std::string> > const& qn)
    {
        assert(qn.size() == 2);
        TwoU1::charge ret;
        ret[0] = details::to_integer (get_quantumnumber(state, qn[0].first) );
        ret[1] = details::to_integer( get_quantumnumber(state, qn[1].first) );
        return ret;
    }

    template <>
    TwoU1::charge init_charge<TwoU1> (const alps::Parameters& parms, std::vector<std::pair<std::size_t, std::string> > const& qn)
    {
        assert(qn.size() == 2);
        TwoU1::charge c = TwoU1::SingletCharge;
        if (parms.defined(qn[0].second+"_total")) {
            alps::half_integer<short> tmp = alps::evaluate<double>(static_cast<std::string>(parms[qn[0].second+"_total"]), parms);
            c[0] = details::to_integer(tmp);
        }
        if (parms.defined(qn[1].second+"_total")) {
            alps::half_integer<short> tmp = alps::evaluate<double>(static_cast<std::string>(parms[qn[1].second+"_total"]), parms);
            c[1] = details::to_integer(tmp);
        }
        return c;
    }

    template <>
    std::map<TwoU1::charge,std::size_t> init_qn_charges<TwoU1>
    (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
    {
        assert(conserved_qn.size() == 2);
        std::map<TwoU1::charge,std::size_t> ret;
        for (int i=0; i<states.size(); ++i)
            ret[convert_alps<TwoU1>(states[i], conserved_qn)] = 1;
        return ret;
    }

    template <>
    std::map<alps::site_state<short>, std::pair<TwoU1::charge, std::size_t> > init_coords<TwoU1>
    (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
    {
        assert(conserved_qn.size() == 2);
        std::map<alps::site_state<short>, std::pair<TwoU1::charge, std::size_t> > ret;
        for (std::size_t i=0; i<states.size(); ++i)
            ret[states[i]] = std::make_pair(convert_alps<TwoU1>(states[i], conserved_qn), 0);
        return ret;
    }

	
}
