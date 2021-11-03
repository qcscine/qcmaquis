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
#include "dmrg/models/generate_mpo/1D_mpo_maker.hpp"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/H2Fixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
> symmetries;

/**
 * @brief Checks that applying the destructor we get the same MPS of the cation.
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE( Test_MPO_Times_MPS_ExpVal, S, symmetries, H2Fixture ) 
{
    // Generates the HF MPS
    parametersH2.set("init_state", "hf");
    parametersH2.set("hf_occ", "4,1");
    auto lattice = Lattice(parametersH2);
    auto model = Model<matrix, S>(lattice, parametersH2);
    auto mpsHF = MPS<matrix, S>(lattice.size(), *(model.initializer(lattice, parametersH2)));
    // Creates the destructor operator
    auto traitClass = MPOTimesMPSTraitClass<matrix, S>(mpsHF, model, lattice, 
                                                       model.total_quantum_numbers(parametersH2),
                                                       parametersH2["max_bond_dimension"]);
    auto ionizedMPS = traitClass.ionizeMPS(0, generate_mpo::IonizedOrbital::Up);
    // Creates the ionized MPS from the mps_intializer
    parametersH2.set("hf_occ", "2,1");
    parametersH2.set("u1_total_charge1", 0);
    parametersH2.set("u1_total_charge2", 1);
    auto modelCation = Model<matrix, S>(lattice, parametersH2);
    auto mpsHFCation = MPS<matrix, S>(lattice.size(), *(modelCation.initializer(lattice, parametersH2)));
    auto overlapBetweenMPS = overlap(mpsHFCation, ionizedMPS);
    // Final check
    BOOST_CHECK_CLOSE(std::abs(overlapBetweenMPS), 1., 1.E-10);
};
