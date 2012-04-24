/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/coded/factory.h"
#include "dmrg/models/coded/lattice.hpp"

inline Lattice_ptr lattice_factory (BaseParameters & parms)
{
    if (parms.get<std::string>("LATTICE") == std::string("periodic chain lattice"))
        return Lattice_ptr(new ChainLattice(parms, true));
    else if (parms.get<std::string>("LATTICE") == std::string("chain lattice"))
        return Lattice_ptr(new ChainLattice(parms, false));
    else if (parms.get<std::string>("LATTICE") == std::string("open chain lattice"))
        return Lattice_ptr(new ChainLattice(parms, false));
    else if (parms.get<std::string>("LATTICE") == std::string("square lattice"))
        return Lattice_ptr(new SquareLattice(parms));
    else if (parms.get<std::string>("LATTICE") == std::string("open square lattice"))
        return Lattice_ptr(new SquareLattice(parms));
    else {
        throw std::runtime_error("Don't know this lattice!");
    }
}
