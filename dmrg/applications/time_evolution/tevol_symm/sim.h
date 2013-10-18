/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/


#include "tevol_sim.h"
#include "tevol_nn_sim.h"
#include "tevol_mpo_sim.h"

template <>
void run_tevol<matrix, grp>(DmrgParameters & parms, ModelParameters & model)
{
    if (parms["te_type"] == "nn")
    {
        tevol_sim<matrix, grp, nearest_neighbors_evolver<matrix, grp> > sim(parms, model);
        sim.run();
    } else if (parms["te_type"] == "mpo") {
        tevol_sim<matrix, grp, mpo_evolver<matrix, grp> > sim(parms, model);
        sim.run();
    } else {
        throw std::runtime_error("Don't know this time evolution. ("+parms["te_type"].str()+")");
    }
}

template <>
void run_tevol<cmatrix, grp>(DmrgParameters & parms, ModelParameters & model)
{
    if (parms["te_type"] == "nn")
    {
        tevol_sim<cmatrix, grp, nearest_neighbors_evolver<cmatrix, grp> > sim(parms, model);
        sim.run();
    } else if (parms["te_type"] == "mpo") {
        tevol_sim<cmatrix, grp, mpo_evolver<cmatrix, grp> > sim(parms, model);
        sim.run();
    } else {
        throw std::runtime_error("Don't know this time evolution. ("+parms["te_type"].str()+")");
    }
}


