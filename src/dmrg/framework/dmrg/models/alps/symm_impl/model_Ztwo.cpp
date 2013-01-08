/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/alps/model_symm.hpp"

	
// Symmetry dependent implementation

// Z2 Symmetry
template <>
Ztwo::charge init_charge<Ztwo> (const alps::Parameters& parms, std::vector<std::pair<std::size_t, std::string> > const& qn)
{
    assert(qn.size() == 1);
    Ztwo::charge c = Ztwo::IdentityCharge;
    if (parms.defined(qn[0].second+"_total")) {
        int tmp = alps::evaluate<double>(static_cast<std::string>(parms[qn[0].second+"_total"]), parms);
        if (!(tmp == 0 || tmp == 1))
            throw std::runtime_error("Invalid value for " + qn[0].second + "_total");
        c = (tmp == 1 ? Ztwo::Minus : Ztwo::Plus);
    }
    return c;
}

template <>
Ztwo::charge convert_alps<Ztwo> (alps::site_state<short> const & state, std::vector<std::pair<std::size_t, std::string> > const& qn)
{
    assert(qn.size() == 1);
    int tmp = detail::to_integer( get_quantumnumber(state, qn[0].first) );
    return (tmp % 2 == 0 ? Ztwo::Plus : Ztwo::Minus);
}

template <>
std::map<Ztwo::charge,std::size_t> init_qn_charges<Ztwo>
(std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
{
    assert(conserved_qn.size() == 1);
    std::map<Ztwo::charge,std::size_t> ret;
    for (int i=0; i<states.size(); ++i)
        ret[convert_alps<Ztwo>(states[i], conserved_qn)] += 1;
    return ret;
}

template <>
std::map<alps::site_state<short>, std::pair<Ztwo::charge, std::size_t> > init_coords<Ztwo>
(std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states)
{
    assert(conserved_qn.size() == 1);
    std::map<alps::site_state<short>, std::pair<Ztwo::charge, std::size_t> > ret;
    std::map<Ztwo::charge, std::size_t> used;
    for (std::size_t i=0; i<states.size(); ++i) {
        Ztwo::charge c = convert_alps<Ztwo>(states[i], conserved_qn);
        ret[states[i]] = std::make_pair(c, used[c]++);
    }
    return ret;
}

