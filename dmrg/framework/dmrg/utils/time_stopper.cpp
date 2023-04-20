/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifdef USE_AMBIENT
#include <ambient/ambient.hpp>
#include <ambient/container/future.hpp>
#endif

#include "dmrg/utils/time_stopper.h"

time_stopper::time_stopper(double timelimit)
: limit(timelimit)
, start(std::chrono::high_resolution_clock::now())
{ }

bool time_stopper::valid() const {
    return limit.count();
}

bool time_stopper::operator()() const {
#ifdef USE_AMBIENT
    ambient::future<double> flag;
    ambient::async( [&](ambient::future<double> & r){
        r.set(limit.count() > 0 && std::chrono::high_resolution_clock::now() > start + limit);
    }, flag );
    return (flag > 0.);
#else
    return (limit.count() > 0 && std::chrono::high_resolution_clock::now() > start + limit);
#endif
}

std::chrono::duration<double> time_stopper::time_left() const {
    return (start + limit) - std::chrono::high_resolution_clock::now();
}
