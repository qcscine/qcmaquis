/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#ifdef DMRG_VIBRATIONAL

#include <random>
#include <boost/filesystem.hpp>
#include <boost/test/included/unit_test.hpp>
#include "Fixtures/EthyleneNModeFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

// Needed for the diagonalization
#include "dmrg/utils/utils.hpp"
#include "utils/timings.h"
#include "utils/traits.hpp"
#include "utils/bindings.hpp"

/**
 * @brief N-mode calculation on two-mode PES of ethylene.
 * The goal of this test is to verify that the SS + noise optimizer works
 * properly. In the original version of the code -- i.e., with the original
 * constructor of the MPS from an ONV -- only the non-zero blocks were added
 * and, as a consequence, the noise was not sufficient to move the optimization
 * from the initial guess. With the new constructor, which also includes the zero
 * blocks of the MPS, the noise is sufficient to move the optimization towards
 * the correct energy. For this reason, we 
 */
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_Ethylene_nMode, EthyleneNModeFixture)
{
#ifdef HAVE_NU1
#if DMRG_NUMSYMM == 12
    parametersEthyleneNMode.set("init_type", "basis_state_generic");
    parametersEthyleneNMode.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
    parametersEthyleneNMode.set("nsweeps", 50);
    parametersEthyleneNMode.set("ngrowsweeps", 10);
    parametersEthyleneNMode.set("nmainsweeps", 10);
    parametersEthyleneNMode.set("max_bond_dimension", 10);
    parametersEthyleneNMode.set("alpha_initial", 1.0E-3);
    parametersEthyleneNMode.set("alpha_main", 1.0E-15);
    parametersEthyleneNMode.set("alpha_final", 0);
    parametersEthyleneNMode.set("truncation_initial", 1.0E-20);
    parametersEthyleneNMode.set("truncation_main", 1.0E-20);
    parametersEthyleneNMode.set("eigensolver", "IETL_JCD");
    parametersEthyleneNMode.set("optimization", "singlesite");
    // Creates the interface
    maquis::DMRGInterface<double> interface(parametersEthyleneNMode);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), 10985.7265979926, 1.0E-5);
#endif // DMRG_NUMSYMM == 12
#endif // HAS_NU1
}

#endif // DMRG_VIBRATIONAL