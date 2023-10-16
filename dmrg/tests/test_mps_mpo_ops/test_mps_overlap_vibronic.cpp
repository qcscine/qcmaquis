/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE MPSOverlapVibronic

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/VibronicFixture.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

#ifdef DMRG_VIBRONIC


/** @brief Checks that overlap between two ONV that are equal but correspond to
 *         different excited states is zero for the exictonic Hamiltonian */
BOOST_FIXTURE_TEST_CASE( Test_MPS_Overlap_Vibronic, VibronicFixture )
{
#ifdef HAVE_U1
    // ONV1 
    parametersExcitonicAggregate.set("vibronic_sorting", "intertwined");
    parametersExcitonicAggregate.set("init_type", "basis_state_generic");
    parametersExcitonicAggregate.set("init_basis_state", "1,1,1,2,3,4,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
    auto excitonicLattice = Lattice(parametersExcitonicAggregate);
    auto excitonicModel = Model<matrix, U1>(excitonicLattice, parametersExcitonicAggregate);
    auto mpsHF1 = MPS<matrix, U1>(excitonicLattice.size(), *(excitonicModel.initializer(excitonicLattice, parametersExcitonicAggregate)));
    // ONV2
    parametersExcitonicAggregate.set("init_type", "basis_state_generic");
    parametersExcitonicAggregate.set("init_basis_state", "0,1,1,2,3,4,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
    excitonicLattice = Lattice(parametersExcitonicAggregate);
    excitonicModel = Model<matrix, U1>(excitonicLattice, parametersExcitonicAggregate);
    auto mpsHF2 = MPS<matrix, U1>(excitonicLattice.size(), *(excitonicModel.initializer(excitonicLattice, parametersExcitonicAggregate)));
    // Overlap calculation
    double overlapHF1 = overlap(mpsHF1, mpsHF2);
    double overlapHF2 = overlap(mpsHF2, mpsHF1);
    BOOST_CHECK_CLOSE(overlapHF1, overlapHF2, 1.E-10);
    BOOST_CHECK_CLOSE(overlapHF1, 0., 1.E-10);
#endif // HAVE_U1
}

#ifdef HAVE_U1

/** @brief Checks Hermitianity of the overlap calculation for the real-valued U1 case */
BOOST_FIXTURE_TEST_CASE( Test_MPS_Overlap_Hermitian_Excitonic, VibronicFixture )
{
    // Global variables
    auto lattice = Lattice(parametersExcitonicAggregate);
    auto model = Model<matrix, U1>(lattice, parametersExcitonicAggregate);
    // Modifies the init parameters to create two different MPSs
    parametersExcitonicAggregate.set("init_type", "const");
    auto mpsConst = MPS<matrix, U1>(lattice.size(), *(model.initializer(lattice, parametersExcitonicAggregate)));
    parametersExcitonicAggregate.set("init_type", "default");
    auto mpsDefault = MPS<matrix, U1>(lattice.size(), *(model.initializer(lattice, parametersExcitonicAggregate)));
    // Calculates the overlap
    double overlapOriginal = overlap(mpsConst, mpsDefault);
    double overlapHerm = overlap(mpsDefault, mpsConst);
    BOOST_CHECK_CLOSE(overlapOriginal, overlapHerm, 1.E-10);
}

#endif // HAVE_U1

#endif // DMRG_VIBRONIC
