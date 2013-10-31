/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "factory_impl.h"

#include "dmrg/models/coded/factory_none.hpp"
#include "dmrg/models/continuum/factory_none.hpp"

typedef TrivialGroup grp;

#if defined USE_AMBIENT
impl_init_model(pmatrix, grp)
impl_init_model(cpmatrix, grp)
#else
impl_init_model(matrix, grp)
impl_init_model(cmatrix, grp)
#endif

template <>
TrivialGroup::charge init_qn<TrivialGroup> (BaseParameters & model)
{
    return TrivialGroup::IdentityCharge;
}
