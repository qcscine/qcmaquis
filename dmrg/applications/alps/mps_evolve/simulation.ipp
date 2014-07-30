/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "dmrg/sim/matrix_types.h"

#include "mps_evolve/tevol_sim.hpp"
#include "dmrg/evolve/tevol_nn_sim.h"
#include "dmrg/evolve/tevol_mpo_sim.h"
#include "mps_evolve/simulation.hpp"

template <class SymmGroup>
void simulation<SymmGroup>::run(DmrgParameters & parms, bool write_xml)
{
    /// Check which matrix to use and which time evolution
    boost::scoped_ptr<abstract_sim> sim;
    if (!parms.defined("COMPLEX") || parms["COMPLEX"]) {
        if (parms["te_type"] == "nn")
            sim.reset(new tevol_sim<cmatrix, SymmGroup, nearest_neighbors_evolver<cmatrix, SymmGroup> >(parms, write_xml));
        else if (parms["te_type"] == "mpo")
            sim.reset(new tevol_sim<cmatrix, SymmGroup, mpo_evolver<cmatrix, SymmGroup> >(parms, write_xml));
    } else {
        if (parms["te_type"] == "nn")
            sim.reset(new tevol_sim<matrix, SymmGroup, nearest_neighbors_evolver<matrix, SymmGroup> >(parms, write_xml));
        else if (parms["te_type"] == "mpo")
            sim.reset(new tevol_sim<matrix, SymmGroup, mpo_evolver<matrix, SymmGroup> >(parms, write_xml));
    }
    
    /// Run
    sim->run();
}
