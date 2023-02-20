/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

/// Adapted from alps/ngs/scheduler/stop_callback.hpp

#ifndef MAQUIS_DMRG__UTIS_TIME_STOPPER_H
#define MAQUIS_DMRG__UTIS_TIME_STOPPER_H

#include <chrono>

class time_stopper {
public:
    time_stopper(double timelimit);
    bool valid() const;
    bool operator()() const;
    std::chrono::duration<double> time_left() const;
private:
    std::chrono::duration<double> limit;
    std::chrono::high_resolution_clock::time_point start;
};

#endif /* defined(MAQUIS_DMRG__UTIS_STOP_CALLBACK_H) */
