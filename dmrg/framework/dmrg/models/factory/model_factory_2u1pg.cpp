/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "model_factory_symm.h"

#include "dmrg/models/coded/factory_2u1pg.hpp"
//#include "dmrg/models/continuum/factory_2u1pg.hpp"

typedef TwoU1PG grp;

#if defined USE_AMBIENT
impl_model_factory(pmatrix, grp)
impl_model_factory(cpmatrix, grp)
#else
impl_model_factory(matrix, grp)
impl_model_factory(cmatrix, grp)
#endif


//template <>
//TwoU1PG::charge init_qn<TwoU1PG> (BaseParameters & model)
//{
//    TwoU1PG::charge initc;
//    initc[0] = model["u1_total_charge1"];
//    initc[1] = model["u1_total_charge2"];
//    initc[2] = model["irrep_charge"];
//    return initc;
//}
