/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_MODELS_H
#define APP_MODELS_H

#include "utils/BaseParameters.h"

#include "lattice.h"
#include "hamiltonian.h"
#include "measurements.h"

namespace app {
    
    template <class Matrix, class SymmGroup>
    void model_parser (std::string const & type, BaseParameters & parms,
                       Lattice* & lattice, Hamiltonian<Matrix, SymmGroup> & H,
                       typename SymmGroup::charge& initc, Measurements<Matrix, SymmGroup>& meas);
}


#endif
