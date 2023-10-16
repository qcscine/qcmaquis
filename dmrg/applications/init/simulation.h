/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_SIM_RUN_H
#define MAQUIS_SIM_RUN_H

#include <boost/shared_ptr.hpp>

#include "dmrg/utils/DmrgParameters.h"

struct simulation_base {
    virtual ~simulation_base() {}
    virtual void run(DmrgParameters & parms) =0;
};

template <class SymmGroup>
struct simulation : public simulation_base {
    void run(DmrgParameters & parms);
};

struct simulation_traits {
    typedef std::shared_ptr<simulation_base> shared_ptr;
    template <class SymmGroup> struct F {
        typedef simulation<SymmGroup> type;
    };
};

#endif
