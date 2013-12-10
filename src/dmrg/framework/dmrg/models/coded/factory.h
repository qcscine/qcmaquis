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

template<class Matrix, class SymmGroup>
struct coded_model_factory {
    static boost::shared_ptr<model_impl<Matrix, SymmGroup> >
    parse(Lattice const &, BaseParameters &);
};

#endif
