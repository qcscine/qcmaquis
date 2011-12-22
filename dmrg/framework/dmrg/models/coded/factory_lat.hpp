/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/coded/factory.h"
#include "dmrg/models/coded/lattice.hpp"

namespace app {
    
    inline Lattice_ptr lattice_factory (BaseParameters & parms)
    {
        if (parms.get<std::string>("LATTICE") == std::string("periodic chain lattice"))
            return Lattice_ptr(new ChainLattice(parms, true));
        else if (parms.get<std::string>("LATTICE") == std::string("open chain lattice"))
            return Lattice_ptr(new ChainLattice(parms, false));
        else {
            throw std::runtime_error("Don't know this lattice!");
        }
    }
    
}
