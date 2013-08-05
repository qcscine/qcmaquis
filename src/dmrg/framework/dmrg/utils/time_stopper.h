/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

/// Adapted from alps/ngs/scheduler/stop_callback.hpp

#ifndef MAQUIS_DMRG__UTIS_TIME_STOPPER_H
#define MAQUIS_DMRG__UTIS_TIME_STOPPER_H

#include <boost/chrono.hpp>

class time_stopper {
public:
    time_stopper(double timelimit);
    bool operator()();
private:
    boost::chrono::duration<double> limit;
    boost::chrono::high_resolution_clock::time_point start;
};

#endif /* defined(MAQUIS_DMRG__UTIS_STOP_CALLBACK_H) */
