/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_CODED_FACTORY_H
#define MAQUIS_DMRG_MODELS_CODED_FACTORY_H

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"

inline Lattice_ptr lattice_factory (BaseParameters &);

template<class Matrix, class SymmGroup>
struct model_factory {
    static typename model_traits<Matrix, SymmGroup>::model_ptr
    parse (Lattice const &, BaseParameters &);
};
    
#include "factory_lat.hpp"
#endif
