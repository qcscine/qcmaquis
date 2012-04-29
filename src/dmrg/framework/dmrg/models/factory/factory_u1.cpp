/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "factory_impl.h"

#include "dmrg/models/coded/factory_u1.hpp"
#include "dmrg/models/continuum/factory_u1.hpp"

typedef U1 grp;

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
U1::charge init_qn<U1> (BaseParameters & model)
{
    return model.get<int>("u1_total_charge");
}
