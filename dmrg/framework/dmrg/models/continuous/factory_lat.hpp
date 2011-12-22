/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/continuous/lattice.hpp"

namespace app {
    
    inline Lattice_ptr cont_lattice_factory (BaseParameters & parms)
    {
        if (parms.get<std::string>("lattice") == std::string("continuous_chain")
            || parms.get<std::string>("lattice") == std::string("continuous_left_chain")
            || parms.get<std::string>("lattice") == std::string("continuous_center_chain"))
            return Lattice_ptr(new ContChain(parms, false));
        else if (parms.get<std::string>("lattice") == std::string("periodic_continuous_chain"))
            return Lattice_ptr(new ContChain(parms, true));
        else {
            throw std::runtime_error("Don't know this lattice!");
        }
    }
    
}
