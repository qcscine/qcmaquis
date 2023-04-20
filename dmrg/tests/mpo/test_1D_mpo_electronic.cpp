/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
    parametersH2.set("init_type", "hf");
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
