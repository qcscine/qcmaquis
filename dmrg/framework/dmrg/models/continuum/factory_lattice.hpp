/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/continuum/lattice.hpp"

inline boost::shared_ptr<lattice_impl> cont_lattice_factory (BaseParameters & parms)
{
    typedef boost::shared_ptr<lattice_impl> impl_ptr;
    std::string lat_str = parms.is_set("lattice") ? "lattice" : "LATTICE";
    if (parms[lat_str] == std::string("continuous_chain")
        || parms[lat_str] == std::string("continuous_left_chain")
        || parms[lat_str] == std::string("continuous_center_chain"))
        return impl_ptr(new ContChain(parms, false));
    else if (parms[lat_str] == std::string("periodic_continuous_chain"))
        return impl_ptr(new ContChain(parms, true));
    else {
        throw std::runtime_error("Don't know this lattice!");
    }
}
