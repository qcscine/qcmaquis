/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/lattice.h"
#include "dmrg/models/coded/factory_lattice.hpp"
#include "dmrg/models/continuum/factory_lattice.hpp"

#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/lattice.hpp"
#endif

#ifdef ENABLE_LL_MODELS
#include "dmrg/models/ll/ll_models.h"
#endif

/// lattice factory
boost::shared_ptr<lattice_impl>
lattice_factory(BaseParameters & parms, BaseParameters & model)
{
    typedef boost::shared_ptr<lattice_impl> impl_ptr;
    
    if (parms["lattice_library"] == "coded") {
        return coded_lattice_factory(parms, model);
    } else if (parms["lattice_library"] == "alps") {
#ifdef ENABLE_ALPS_MODELS
        return impl_ptr( new alps_lattice(model) );
#else
        throw std::runtime_error("This code was compiled without alps lattice.");
#endif
    } else if (parms["lattice_library"] == "continuum") {
        return cont_lattice_factory(model);
#ifdef ENABLE_LL_MODELS
    } else if (parms["lattice_library"] == "ll") {
        return ll_lattice_factory(model);
#endif
    } else {
        throw std::runtime_error("Don't know this lattice_library!");
    }
}

