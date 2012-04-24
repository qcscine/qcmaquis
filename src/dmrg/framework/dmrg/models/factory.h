/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_MODELS_H
#define APP_MODELS_H

#include "dmrg/utils/BaseParameters.h"

#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"

template <class Matrix, class SymmGroup>
void model_parser (std::string lattice_lib, std::string model_lib,
                   BaseParameters & parms,
                   Lattice_ptr & lattice,
                   typename model_traits<Matrix, SymmGroup>::model_ptr & model);

#endif
