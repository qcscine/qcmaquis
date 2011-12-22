/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_CONT_FACTORY_H
#define MAQUIS_DMRG_MODELS_CONT_FACTORY_H

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"

namespace app {
    
    inline Lattice_ptr cont_lattice_factory (BaseParameters &);
    
    
    template<class Matrix, class SymmGroup>
    struct cont_model_factory {
        static typename model_traits<Matrix, SymmGroup>::model_ptr
        parse (Lattice const &, BaseParameters &);
    };
    
}

#include "factory_lat.hpp"

#endif
