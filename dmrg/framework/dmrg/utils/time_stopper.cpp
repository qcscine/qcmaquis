/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/utils/time_stopper.h"

time_stopper::time_stopper(double timelimit)
: limit(timelimit)
, start(boost::chrono::high_resolution_clock::now())
{ }

bool time_stopper::operator()() {
    return (limit.count() > 0 && boost::chrono::high_resolution_clock::now() > start + limit);
}
