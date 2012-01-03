/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/


#include "tevol_sim.h"

template <>
void run_tevol<matrix, grp>(DmrgParameters & parms, ModelParameters & model)
{
    dmrg_tevol_sim<matrix, grp> sim(parms, model);
    sim.run();
}

template <>
void run_tevol<cmatrix, grp>(DmrgParameters & parms, ModelParameters & model)
{
    dmrg_tevol_sim<cmatrix, grp> sim(parms, model);
    sim.run();
}


