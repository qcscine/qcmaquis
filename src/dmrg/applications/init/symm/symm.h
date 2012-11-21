/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/


#include "sim.h"

namespace maquis { namespace dmrg {
    template <>
    void run_sim<grp>(DmrgParameters & parms, ModelParameters & model)
    {
        dmrg_init<matrix, grp> sim(parms, model);
        sim.build();
    }
} }
