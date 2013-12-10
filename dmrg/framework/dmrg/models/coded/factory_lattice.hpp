/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_CODED_FACTORY_LATTICE_H
#define MAQUIS_DMRG_MODELS_CODED_FACTORY_LATTICE_H

#include "dmrg/models/coded/lattice.hpp"

inline boost::shared_ptr<lattice_impl>
coded_lattice_factory(BaseParameters & parms, BaseParameters & model)
{
    typedef boost::shared_ptr<lattice_impl> impl_ptr;
    if (model["LATTICE"] == std::string("periodic chain lattice"))
        return impl_ptr(new ChainLattice(model, true));
    else if (model["LATTICE"] == std::string("chain lattice"))
        return impl_ptr(new ChainLattice(model, false));
    else if (model["LATTICE"] == std::string("open chain lattice"))
        return impl_ptr(new ChainLattice(model, false));
    else if (model["LATTICE"] == std::string("square lattice"))
        return impl_ptr(new SquareLattice(model));
    else if (model["LATTICE"] == std::string("open square lattice"))
        return impl_ptr(new SquareLattice(model));
    else if (model["LATTICE"] == std::string("orbitals"))
        return impl_ptr(new Orbitals(model));
    else {
        throw std::runtime_error("Don't know this lattice!");
    }
}

#endif
