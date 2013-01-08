/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "factory_impl.h"

#include "dmrg/models/coded/factory_Ztwo.hpp"
#include "dmrg/models/continuum/factory_Ztwo.hpp"

typedef Ztwo grp;

#if defined USE_MTM
impl_init_model(mtmatrix, grp)
impl_init_model(cmtmatrix, grp)
#elif defined USE_AMBIENT
impl_init_model(pmatrix, grp)
impl_init_model(cpmatrix, grp)
#else
impl_init_model(matrix, grp)
impl_init_model(cmatrix, grp)
#endif

template <>
Ztwo::charge init_qn<Ztwo> (BaseParameters & model)
{
    return Ztwo::IdentityCharge;
}
