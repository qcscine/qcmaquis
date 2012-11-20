/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/


#include "dmrg_sim.h"

namespace maquis { namespace dmrg {
    template <>
    void run_sim<grp>(DmrgParameters & parms, ModelParameters & model)
    {
        dmrg_sim<matrix, grp> sim(parms, model);
        sim.measure();
    }
} }
