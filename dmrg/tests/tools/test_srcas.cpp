/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE SRCAS

#include "dmrg/tools/srcas_utilities.h"

#include "Fixtures/NModeFixture.h"
#include "Fixtures/WatsonFixture.h"
#include "Fixtures/LiHFixture.h"

#include <numeric>
#include <algorithm>
#include <memory>
#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/list.hpp>

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_TwoU1
, TwoU1
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
#ifdef HAVE_SU2U1
, SU2U1
#endif
> symmetries;

/** @brief Test SRCAS for electronic calculations for all possible symmetries */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_LiH_DMRG_SSvsTS, S, symmetries, LiHFixture)
{
    using InterfaceType = maquis::DMRGInterface<double>;
    // Generic parameters
    parametersLiH.set("max_bond_dimension", 50);
    parametersLiH.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
    parametersLiH.set("init_type", "const");
    parametersLiH.set("nsweeps", 20);
    parametersLiH.set("ngrowsweeps", 2);
    parametersLiH.set("nmainsweeps", 5);
    parametersLiH.set("optimizer", "singlesite");
    parametersLiH.set("hf_occ", "4,1,1,1");
    parametersLiH.set("srcas_numSamples", 1000);
    parametersLiH.set("srcas_overlapThreshold", -1);
    parametersLiH.set("spin", 0);
    // Creates the interface and performs a optimization so that we have a MPS to compare to
    std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(parametersLiH);
    interface->optimize();
    // Creates SRCAS object
    SRCAS<double> srcas(parametersLiH, interface);
    std::vector<int> currQueen = srcas.getCurrentQueen();
    BOOST_CHECK_EQUAL(currQueen.size(), parametersLiH["L"]);
    int sumOfQueen = std::accumulate(currQueen.begin(), currQueen.end(), 0);
    BOOST_CHECK_EQUAL(sumOfQueen, 7);
    srcas.run();
    std::map<std::vector<int>, double> detTable = srcas.getDetTable();
    srcas.printResults();
    // Check that all dets have been sampled
    if (parametersLiH["symmetry"]=="2u1pg" || parametersLiH["symmetry"]=="2u1")
        BOOST_CHECK_EQUAL(detTable.size(), 16);
    else if (parametersLiH["symmetry"]=="su2u1pg" || parametersLiH["symmetry"]=="su2u1")
        BOOST_CHECK_EQUAL(detTable.size(), 28);

    BOOST_CHECK_CLOSE(srcas.getCompleteness(), 1.0, 1.0E-10); // Can be this thight, because we should sample all dets

    // Perform an additional test for the spin symmetry for the triplet state
    if (parametersLiH["symmetry"]=="su2u1pg" || parametersLiH["symmetry"]=="su2u1") {
        parametersLiH.set("spin", 2);
        // Creates the interface and performs a optimization so that we have a MPS to compare to
        std::shared_ptr<InterfaceType> interface2 = std::make_shared<InterfaceType>(parametersLiH);
        interface2->optimize();
        // Creates SRCAS object
        SRCAS<double> srcas2(parametersLiH, interface2);
        std::vector<int> currQueen2 = srcas2.getCurrentQueen();
        BOOST_CHECK_EQUAL(currQueen2.size(), parametersLiH["L"]);
        srcas2.run();
        std::map<std::vector<int>, double> detTable2 = srcas2.getDetTable();
        srcas2.printResults();
        BOOST_CHECK_EQUAL(detTable2.size(), 24);
        BOOST_CHECK_CLOSE(srcas2.getCompleteness(), 1.0, 1.0E-10); // Can be this thight, because we should sample all dets
    }
}

#ifdef HAVE_TrivialGroup

/**
 * @brief SRCAS test for watson-based harmonic ethylene calculation
 */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_SRCAS_Ethylene_Harmonic, WatsonFixture)
{
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersEthyleneWatsonHarmonic.set("init_type", "basis_state_generic");
    parametersEthyleneWatsonHarmonic.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
    parametersEthyleneWatsonHarmonic.set("nsweeps", 5);
    parametersEthyleneWatsonHarmonic.set("max_bond_dimension", 20);
    parametersEthyleneWatsonHarmonic.set("MODEL", "watson");
    // Creates the interface and performs a optimization so that we have a MPS to compare to
    std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(parametersEthyleneWatsonHarmonic);
    interface->optimize();
    // Creates SRCAS object
    SRCAS<double> srcas(parametersEthyleneWatsonHarmonic, interface);
    std::vector<int> currQueen = srcas.getCurrentQueen();
    BOOST_CHECK_EQUAL(currQueen.size(), parametersEthyleneWatsonHarmonic["L"]);
    int sumOfQueen = std::accumulate(currQueen.begin(), currQueen.end(), 0);
    BOOST_CHECK_EQUAL(sumOfQueen, 0);
    srcas.run();
    std::map<std::vector<int>, double> detTable = srcas.getDetTable();
    BOOST_CHECK_EQUAL(detTable.size(), 1);
    BOOST_CHECK_CLOSE(srcas.getCompleteness(), 1.0, 1.0E-10); // Can be this thight, because we start in GS and have no coupling
}

/**
 * @brief SRCAS test based on Watson-based calculation on the ethylene PES
 * The PES has been taken from the vHBCI reference work, which is
 * J. Chem. Phys., 154, 074104 (2021).
 * Note that the calculation uses the single-site optimizer.
 */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_SRCAS_Ethylene_Sextic_SingleSite, WatsonFixture)
{
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersEthyleneWatson.set("init_type", "basis_state_generic");
    parametersEthyleneWatson.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
    parametersEthyleneWatson.set("nsweeps", 20);
    parametersEthyleneWatson.set("max_bond_dimension", 50);
    parametersEthyleneWatson.set("MODEL", "watson");
    parametersEthyleneWatson.set("optimization", "singlesite");
    parametersEthyleneWatson.set("ngrowsweeps", 2);
    parametersEthyleneWatson.set("nmainsweeps", 2);
    parametersEthyleneWatson.set("alpha_initial", 1.0E-8);
    parametersEthyleneWatson.set("alpha_main", 1.0E-15);
    parametersEthyleneWatson.set("alpha_final", 0.);
    // Creates the interface
    std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(parametersEthyleneWatson);
    interface->optimize();
    // Creates SRCAS object
    SRCAS<double> srcas(parametersEthyleneWatson, interface);
    std::vector<int> currQueen = srcas.getCurrentQueen();
    BOOST_CHECK_EQUAL(currQueen.size(), parametersEthyleneWatson["L"]);
    int sumOfQueen = std::accumulate(currQueen.begin(), currQueen.end(), 0);
    BOOST_CHECK_EQUAL(sumOfQueen, 0);
    srcas.run();
    // Check again now that we are on a different queen
    currQueen = srcas.getCurrentQueen();
    BOOST_CHECK_EQUAL(currQueen.size(), parametersEthyleneWatson["L"]);
    std::map<std::vector<int>, double> detTable = srcas.getDetTable();
    BOOST_CHECK_EQUAL(detTable.size(), 140);
    BOOST_CHECK_CLOSE(srcas.getCompleteness(), 1.0, 1.0E-01);
    srcas.printResults();
}

#endif    

#ifdef HAVE_NU1

/**
 * @brief SRCAS test based on N-mode calculation on one-mode PES of FAD.
 * The basis set has been generated with a DVR primitive basis.
 * Note that we also check the excitation energy.
 */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_SRCAS_FAD_1ModeHamiltonian_ExcitedState, NModeFixture)
{
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersFADOneBodyBinary.set("init_type", "basis_state_generic");
    parametersFADOneBodyBinary.set("init_basis_state", "0");
    parametersFADOneBodyBinary.set("nsweeps", 2);
    parametersFADOneBodyBinary.set("max_bond_dimension",100);
    parametersFADOneBodyBinary.set("chkpfile", "GS.checkpoint.h5");
    parametersFADOneBodyBinary.set("resfule", "GS.results.h5");
    // Creates the interface and checks the resulting energy
    std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(parametersFADOneBodyBinary);
    interface->optimize();
    // Checks that the overlap of the final wave function with the hf determinant is = 1.
    auto targetOverlap = interface->getCICoefficient("0");
    BOOST_CHECK_CLOSE(std::abs(targetOverlap), 1.0, 1.0E-16);
    // Checks that overlap calculations cannot be performed for invalid reference determinants
    BOOST_CHECK_THROW(interface->getCICoefficient("0,0,0"), std::runtime_error);
    // Creates SRCAS object
    SRCAS<double> srcas(parametersFADOneBodyBinary, interface);
    std::vector<int> currQueen = srcas.getCurrentQueen();
    BOOST_CHECK_EQUAL(currQueen.size(), parametersFADOneBodyBinary["nmode_num_modes"]);
    int sumOfQueen = std::accumulate(currQueen.begin(), currQueen.end(), 0);
    BOOST_CHECK_EQUAL(sumOfQueen, 0);
    srcas.run();
    std::map<std::vector<int>, double> detTable = srcas.getDetTable();
    BOOST_CHECK_EQUAL(detTable.size(), 1);
    BOOST_CHECK_CLOSE(srcas.getCompleteness(), 1.00, 1.0E-10); // Can be this thight, because we start in GS and have no coupling
    // Excited-state calculations
    parametersFADOneBodyBinary.set("chkpfile", "ES.checkpoint.h5");
    parametersFADOneBodyBinary.set("resfule", "ES.results.h5");
    parametersFADOneBodyBinary.set("nsweeps", 20);
    parametersFADOneBodyBinary.set("init_type", "const");
    parametersFADOneBodyBinary.set("n_ortho_states", 1);
    parametersFADOneBodyBinary.set("ortho_states", "GS.checkpoint.h5");
    // Creates a new interface object and reruns the optimization
    std::shared_ptr<InterfaceType> interfaceES = std::make_shared<InterfaceType>(parametersFADOneBodyBinary);
    interfaceES->optimize();
    auto esEnergy = interfaceES->energy();
    // Creates SRCAS object for excited state
    parametersFADOneBodyBinary.set("init_type", "basis_state_generic");
    parametersFADOneBodyBinary.set("init_basis_state", "1");
    SRCAS<double> srcasES(parametersFADOneBodyBinary, interfaceES);
    currQueen = srcasES.getCurrentQueen();
    BOOST_CHECK_EQUAL(currQueen.size(), parametersFADOneBodyBinary["nmode_num_modes"]);
    sumOfQueen = std::accumulate(currQueen.begin(), currQueen.end(), 0);
    BOOST_CHECK_EQUAL(sumOfQueen, 1);
    srcasES.run();
    std::map<std::vector<int>, double> detTableES = srcasES.getDetTable();
    auto maxAbsCICoeff = std::max_element(detTableES.begin(), detTableES.end(), []
        (const std::pair<std::vector<int>, double> a, const std::pair<std::vector<int>, double> b)
        {return std::abs(a.second) < std::abs(b.second);});
    BOOST_CHECK_CLOSE(std::abs(maxAbsCICoeff->second), 1.0, 1.0E-3);
    BOOST_CHECK_EQUAL(std::accumulate(maxAbsCICoeff->first.begin(), maxAbsCICoeff->first.end(), 0), 1);
    BOOST_CHECK_CLOSE(srcasES.getCompleteness(), 1.00, 1.0E-5);
    boost::filesystem::remove_all("GS.results.h5");
    boost::filesystem::remove_all("GS.checkpoint.h5");
    boost::filesystem::remove_all("ES.results.h5");
    boost::filesystem::remove_all("ES.checkpoint.h5");
}

#endif
