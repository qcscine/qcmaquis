/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/


#include "mg_meas_sim.h"

template <>
void run_mg_meas<grp>(DmrgParameters & parms, ModelParameters & model)
{
    mg_meas_sim<Matrix, grp> sim(parms, model);
    sim.run();
}

