/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/sim/matrix_types.h"

#include "mg_meas_sim.h"
#include "mg_meas_traits.h"

template <class SymmGroup>
void simulation<SymmGroup>::run(DmrgParameters & parms)
{
    mg_meas_sim<matrix, SymmGroup> sim(parms);
    sim.run();
}
