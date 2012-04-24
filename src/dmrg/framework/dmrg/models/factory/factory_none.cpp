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

#ifdef USE_MTM
impl_init_model(mtmatrix1, grp)
impl_init_model(cmtmatrix1, grp)
#else
impl_init_model(matrix1, grp)
impl_init_model(cmatrix1, grp)
#endif

template <>
TrivialGroup::charge init_qn<TrivialGroup> (BaseParameters & model)
{
    return TrivialGroup::IdentityCharge;
}
