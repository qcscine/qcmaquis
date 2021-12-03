/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Robin Feldmann <robin.feldmann@phys.chem.ethz.ch>
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

#define BOOST_TEST_MODULE MODEL_VIBRONIC_U1

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/vibrational/u1/VibronicModel.hpp"
#include "Fixtures/VibronicFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

/** Test for the integral parser for the trivial vibronic Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_Vibronic, VibronicFixture)
{
#ifdef HAVE_U1
    auto lattice = Lattice(parametersFakeVibronic);
    auto integrals = Vibrational::detail::parseIntegralVibronic<double>(parametersFakeVibronic, lattice);
    // Checks sizes
    BOOST_CHECK_EQUAL(integrals.first.size(), 6);
    BOOST_CHECK_EQUAL(integrals.second.size(), 6);
    BOOST_CHECK_EQUAL(integrals.first[0].size(), 4);
#endif // HAVE_U1
}

#ifdef HAVE_U1


/** Tests the [create_terms] method for the fake vibronic model */
BOOST_FIXTURE_TEST_CASE(Test_VibronicModel_Create_Terms, VibronicFixture)
{
    auto lattice = Lattice(parametersFakeVibronic);
    auto fakeVibronicModel = VibronicModel<tmatrix<double>>(lattice, parametersFakeVibronic);
    auto sizeBefore = fakeVibronicModel.hamiltonian_terms().size();
    BOOST_CHECK_EQUAL(sizeBefore, 0);
    fakeVibronicModel.create_terms();
    auto sizeAfter = fakeVibronicModel.hamiltonian_terms().size();
    BOOST_CHECK_EQUAL(sizeAfter, 6);
}

#endif // HAVE_U1

///** Checks consistency for the physical dimensions for a 1-mode system */
//BOOST_FIXTURE_TEST_CASE(Test_Model_PhysDim_OneMode, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADOneBody);
//    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
//    const auto& physicalDimensions0 = nModeModel.phys_dim(0);
//    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), 2);
//}
//
///** Checks consistency for the physical dimensions for a 2-mode system */
//BOOST_FIXTURE_TEST_CASE(Test_Model_PhysDim_TwoMode, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADTwoBody);
//    auto nModeModel = NMode<tmatrix<double>, 2>(lattice, parametersFADTwoBody, false);
//    const auto& physicalDimensions0 = nModeModel.phys_dim(0);
//    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), 2);
//    const auto& physicalDimensions1 = nModeModel.phys_dim(1);
//    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), 2);
//}
//
///** Checks consistency for the overall QN for a 1-mode system */
//BOOST_FIXTURE_TEST_CASE(Test_Model_TotalQN_OneMode, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADOneBody);
//    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
//    auto totalQN = nModeModel.total_quantum_numbers(parametersFADOneBody);
//    BOOST_CHECK_EQUAL(totalQN[0], 1);
//}
//
///** Checks consistency for the overall QN for a 1-mode system */
//BOOST_FIXTURE_TEST_CASE(Test_Model_TotalQN_TwoMode, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADTwoBody);
//    auto nModeModel = NMode<tmatrix<double>, 2>(lattice, parametersFADTwoBody, false);
//    auto totalQN = nModeModel.total_quantum_numbers(parametersFADTwoBody);
//    BOOST_CHECK_EQUAL(totalQN[0], 1);
//    BOOST_CHECK_EQUAL(totalQN[1], 1);
//    // Now we use a larger lattice
//    auto nModeModelWrongDim = NMode<tmatrix<double>, 3>(lattice, parametersFADTwoBody, false);
//    auto totalQNWrongDim = nModeModelWrongDim.total_quantum_numbers(parametersFADTwoBody);
//    BOOST_CHECK_EQUAL(totalQNWrongDim[0], 1);
//    BOOST_CHECK_EQUAL(totalQNWrongDim[1], 1);
//    BOOST_CHECK_EQUAL(totalQNWrongDim[2], 0);
//}
//
///** Simple check on tags */
//BOOST_FIXTURE_TEST_CASE(Test_Model_Tag_SimpleCheck_OneMode, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADOneBody);
//    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
//    auto identityTag = nModeModel.identity_matrix_tag(0);
//    auto fillingTag = nModeModel.filling_matrix_tag(0);
//    // The nMode Hamiltonian is bosonic, so the tag should be the same
//    BOOST_CHECK(identityTag == fillingTag);
//}
//
///** Simple check on tags for the two-mode Hamiltonian */
//BOOST_FIXTURE_TEST_CASE(Test_Model_Tag_SimpleCheck_TwoMode, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADTwoBody);
//    auto nModeModel = NMode<tmatrix<double>, 2>(lattice, parametersFADTwoBody, false);
//    auto identityTag = nModeModel.filling_matrix_tag(0);
//    auto fillingTag = nModeModel.filling_matrix_tag(1);
//    BOOST_CHECK(identityTag != fillingTag);
//}
//
///** Check on symbolic operator getter */
//BOOST_FIXTURE_TEST_CASE(Test_Model_Symbolic_Operator_OneMode, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADOneBody);
//    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
//    BOOST_CHECK(nModeModel.filling_matrix_tag(0) == nModeModel.get_operator_tag("fill", 0));
//}
//
///** Check on symbolic operator getter for a two-mode Hamiltonian */
//BOOST_FIXTURE_TEST_CASE(Test_Model_Symbolic_Operator_TwoMode, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADOneBody);
//    auto nModeModel = NMode<tmatrix<double>, 2>(lattice, parametersFADOneBody, false);
//    BOOST_CHECK(nModeModel.filling_matrix_tag(0) == nModeModel.get_operator_tag("fill", 0));
//    BOOST_CHECK(nModeModel.identity_matrix_tag(1) == nModeModel.get_operator_tag("id", 1));
//}
//
///** Check on operator table getter */
//BOOST_FIXTURE_TEST_CASE(Test_Model_Operator_Table, NModeFixture)
//{
//    auto lattice = Lattice(parametersFADOneBody);
//    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
//    auto createTag = nModeModel.get_operator_tag("bdag", 0);
//    auto destroyTag = nModeModel.get_operator_tag("b", 0);
//    auto countTag = nModeModel.get_operator_tag("n", 0);
//    auto countTagFromTable = nModeModel.operators_table()->get_product_tag(destroyTag, createTag);
//    auto siteOperatorCount = nModeModel.operators_table()->get_op(countTag);
//    auto siteOperatorCountFromTable = nModeModel.operators_table()->get_op(countTagFromTable.first);
//    BOOST_CHECK_CLOSE(siteOperatorCount.norm(), 1., 1.0E-16);
//    BOOST_CHECK_CLOSE(siteOperatorCountFromTable.norm(), 1., 1.0E-16);
//    auto differenceOperator = siteOperatorCount - siteOperatorCountFromTable;
//    BOOST_CHECK_CLOSE(differenceOperator.norm(), 0., 1.0E-16);
//}
//
//
//#endif // HAVE_NU1
