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
TrivialGroup::charge init_qn<TrivialGroup> (BaseParameters & model)
{
    return TrivialGroup::IdentityCharge;
}
