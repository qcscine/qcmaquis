/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/


#include "dmrg_sim.h"

template <>
void run_dmrg<grp>(DmrgParameters & parms, ModelParameters & model)
{
    dmrg_sim<Matrix, grp> sim(parms, model);
    sim.measure();
}

