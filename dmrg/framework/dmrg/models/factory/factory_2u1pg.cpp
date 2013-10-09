/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "factory_impl.h"

#include "dmrg/models/coded/factory_2u1pg.hpp"
#include "dmrg/models/continuum/factory_2u1pg.hpp"

typedef TwoU1PG grp;

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
TwoU1PG::charge init_qn<TwoU1PG> (BaseParameters & model)
{
    TwoU1PG::charge initc;
    initc[0] = model["u1_total_charge1"];
    initc[1] = model["u1_total_charge2"];
    return initc;
}
