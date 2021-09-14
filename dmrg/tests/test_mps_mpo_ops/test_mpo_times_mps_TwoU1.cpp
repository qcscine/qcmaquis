/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021 Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/BenzeneFixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
> symmetries;

/**
 * @brief Checks that [mpo_times_mps] gives results that are coherent with expval.
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE( Test_MPO_Times_MPS_ExpVal, S, symmetries, BenzeneFixture ) 
{
    // Generates the HF MPS
    auto lattice = Lattice(parametersBenzene);
    auto modelHF = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsHF = MPS<matrix, S>(lattice.size(), *(modelHF.initializer(lattice, parametersBenzene)));
    auto mpo = make_mpo(lattice, modelHF);
    // Calculates the MPS-MPO contraction
    int max_site_type = 0;
    std::vector<int> site_types(lattice.size(), 0);
    for (int p = 0; p < lattice.size(); ++p) {
      site_types[p] = lattice.template get_prop<int>("type", p);
      max_site_type = std::max(site_types[p], max_site_type);
    }
    std::vector<Index<S> > site_bases(max_site_type+1);
    for (int type = 0; type < site_bases.size(); ++type)
      site_bases[type] = modelHF.phys_dim(type);
    auto totalQN = modelHF.total_quantum_numbers(parametersBenzene);
    std::vector<typename S::charge> charges = {S::IdentityCharge};
    std::map<int, Index<S> > mapTrackingBlocks;
    Index<S> tmp;
    tmp.insert(std::make_pair(S::IdentityCharge, 1));
    mapTrackingBlocks[0] = tmp;
    MPS<matrix, S> outputMPS(lattice.size());
    auto indexAllowed = allowed_sectors(site_types, site_bases, totalQN, parametersBenzene["max_bond_dimension"]);
    for (int iMPS = 0; iMPS < mpsHF.length(); iMPS++)
      outputMPS[iMPS] = MPOTimesMPSTraitClass<matrix, matrix, S>::mpo_times_mps(mpo, mpsHF, iMPS, charges, mapTrackingBlocks, indexAllowed);
    // Calculates the energy in two ways
    auto energyFromMPSTimesMPO = overlap(mpsHF, outputMPS)/norm(mpsHF) + mpo.getCoreEnergy();
    auto energyFromExpVal = expval(mpsHF, mpo)/norm(mpsHF);
    BOOST_CHECK_CLOSE(energyFromMPSTimesMPO, energyFromExpVal, 1.E-10);
};