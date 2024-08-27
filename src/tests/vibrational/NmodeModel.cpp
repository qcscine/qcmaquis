/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE MODEL_VIBRATIONAL_NONE

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/vibrational/none/NModeModelCompact.hpp"
#include "dmrg/models/vibrational/none/NModeModelPairedOperators.hpp"
#include "dmrg/models/vibrational/nu1/model.hpp"
#include "Fixtures/NModeFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

BOOST_FIXTURE_TEST_CASE(Test_NmodeModels_Energy, NModeFixture)
{
#ifdef DMRG_VIBRATIONAL
#ifdef HAVE_TrivialGroup
    parametersWater.set("nmode_num_basis", "7,7,7");
    parametersWater.set("init_basis_state", "1,1,1");
    parametersWater.set("symmetry", "none");
    parametersWater.set("integral_file", "integral_file_test_OneBodyWater");
    parametersWater.set("MODEL", "nmodecompact");
    parametersWater.set("LATTICE", "watson lattice");
    parametersWater.set("L", 3);
    maquis::DMRGInterface<std::complex<double>> interface_compact_1body(parametersWater);
    parametersWater.set("nmode_num_basis", "6,6,6");
    parametersWater.set("init_basis_state", "0,0,0");
    parametersWater.set("symmetry", "none");
    parametersWater.set("integral_file", "integral_file_test_OneBodyWater");
    parametersWater.set("MODEL", "nmodecompactpaired");
    parametersWater.set("LATTICE", "watson lattice");
    parametersWater.set("L", 3);
    maquis::DMRGInterface<std::complex<double>> interface_paired_1body(parametersWater);
    parametersWater.set("nmode_num_basis", "6,6,6");
    parametersWater.set("init_basis_state", "0,0,0");
    parametersWater.set("symmetry", "nu1");
    parametersWater.set("integral_file", "integral_file_test_OneBodyWater");
    parametersWater.set("MODEL", "nmode");
    parametersWater.set("LATTICE", "nmode lattice");
    parametersWater.set("L", 18);
    maquis::DMRGInterface<std::complex<double>> interface_nu1_1body(parametersWater);
    interface_compact_1body.optimize();
    interface_paired_1body.optimize();
    interface_nu1_1body.optimize();
    double interface_compact_optimizedEnergy_1body = interface_compact_1body.energy().real();
    double interface_paired_optimizedEnergy_1body = interface_paired_1body.energy().real();
    double interface_nu1_optimizedEnergy_1body = interface_nu1_1body.energy().real();
    BOOST_CHECK_CLOSE(interface_compact_optimizedEnergy_1body, interface_paired_optimizedEnergy_1body, 1.0E-8);
    BOOST_CHECK_CLOSE(interface_compact_optimizedEnergy_1body, interface_nu1_optimizedEnergy_1body, 1.0E-8);
    BOOST_CHECK_CLOSE(interface_paired_optimizedEnergy_1body, interface_nu1_optimizedEnergy_1body, 1.0E-8);
    BOOST_CHECK_CLOSE(interface_compact_optimizedEnergy_1body, 4725.63505369 , 1.0E-8);
    BOOST_CHECK_CLOSE(interface_paired_optimizedEnergy_1body, 4725.63505369, 1.0E-8);
    BOOST_CHECK_CLOSE(interface_nu1_optimizedEnergy_1body, 4725.63505369, 1.0E-8);
#endif //DMRG_VIBRATIONAL
#endif // HAVE_TrivialGroup
}