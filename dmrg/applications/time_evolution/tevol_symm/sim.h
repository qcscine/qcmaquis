/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/


#include "tevol_nn_sim.h"
#include "tevol_mpo_sim.h"

template <>
void run_tevol<matrix, grp>(DmrgParameters & parms, ModelParameters & model)
{
    if (parms["te_type"] == "nn")
    {
        dmrg_tevol_nn_sim<matrix, grp> sim(parms, model);
        sim.run();
    } else if (parms["te_type"] == "mpo") {
        dmrg_tevol_mpo_sim<matrix, grp> sim(parms, model);
        sim.run();
    } else {
        throw std::runtime_error("Don't know this time evolution. ("+parms["te_type"]+")");
    }
}

template <>
void run_tevol<cmatrix, grp>(DmrgParameters & parms, ModelParameters & model)
{
    if (parms["te_type"] == "nn")
    {
        dmrg_tevol_nn_sim<cmatrix, grp> sim(parms, model);
        sim.run();
    } else if (parms["te_type"] == "mpo") {
        dmrg_tevol_mpo_sim<cmatrix, grp> sim(parms, model);
        sim.run();
    } else {
        throw std::runtime_error("Don't know this time evolution. ("+parms["te_type"]+")");
    }
}


