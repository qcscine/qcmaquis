/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "model_factory_symm.h"

#include "dmrg/models/coded/factory_2u1.hpp"
//#include "dmrg/models/continuum/factory_2u1.hpp"

typedef TwoU1 grp;

#if defined USE_AMBIENT
impl_model_factory(pmatrix, grp)
impl_model_factory(cpmatrix, grp)
#else
impl_model_factory(matrix, grp)
impl_model_factory(cmatrix, grp)
#endif
