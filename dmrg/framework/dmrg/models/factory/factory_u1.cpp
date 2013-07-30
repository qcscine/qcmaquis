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

#if defined USE_AMBIENT
template struct cont_model_factory<pmatrix, grp>;
template struct cont_model_factory<cpmatrix, grp>;
impl_init_model(pmatrix, grp)
impl_init_model(cpmatrix, grp)
#else
template struct cont_model_factory<matrix, grp>;
template struct cont_model_factory<cmatrix, grp>;
impl_init_model(matrix, grp)
impl_init_model(cmatrix, grp)
#endif

template <>
U1::charge init_qn<U1> (BaseParameters & model)
{
    return model["u1_total_charge"];
}
