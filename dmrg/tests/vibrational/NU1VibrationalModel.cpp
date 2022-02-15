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

#define BOOST_TEST_MODULE MODEL_VIBRATIONAL_NU1

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/vibrational/nu1/model.hpp"
#include "Fixtures/NModeFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

/** Test for the integral parser with the one-body Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_OneBody, NModeFixture)
{
#ifdef HAVE_NU1
    auto integrals = Vibrational::detail::NModeIntegralParser<double>(parametersFADOneBody, lattice_factory(parametersFADOneBody));
    // Checks sizes
    BOOST_CHECK_EQUAL(integrals.first.size(), 78);
    BOOST_CHECK_EQUAL(integrals.second.size(), 78);
    BOOST_CHECK_EQUAL(integrals.first[0].size(), 12);
    // Check numeric values
    BOOST_CHECK_CLOSE(-2.359242429009664e+03, integrals.second[0], 1e-12);
    BOOST_CHECK_EQUAL(1, integrals.first[0][0]);
    BOOST_CHECK_EQUAL(0, integrals.first[0][1]);
    BOOST_CHECK_EQUAL(1, integrals.first[0][2]);
    BOOST_CHECK_EQUAL(0, integrals.first[0][3]);
    BOOST_CHECK_EQUAL(-1, integrals.first[0][4]);
    BOOST_CHECK_EQUAL(-1, integrals.first[0][5]);
    BOOST_CHECK_EQUAL(-1, integrals.first[0][6]);
    BOOST_CHECK_EQUAL(-1, integrals.first[0][7]);
    BOOST_CHECK_EQUAL(-1, integrals.first[0][8]);
    BOOST_CHECK_EQUAL(-1, integrals.first[0][9]);
    BOOST_CHECK_EQUAL(-1, integrals.first[0][10]);
    BOOST_CHECK_EQUAL(-1, integrals.first[0][11]);
#endif // HAVE_NU1
}

#ifdef HAVE_NU1

/** Test for the integral parser with the two-body Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_TwoBody, NModeFixture)
{
    auto integrals = Vibrational::detail::NModeIntegralParser<double>(parametersFADTwoBody, lattice_factory(parametersFADTwoBody));
    // Checks sizes
    BOOST_CHECK_EQUAL(integrals.first.size(), 4845);
    BOOST_CHECK_EQUAL(integrals.second.size(), 4845);
    // Check numeric values
    int idx = 4844;
    BOOST_CHECK_CLOSE(2.904591355839877e+04, integrals.second[idx], 1e-12);
    BOOST_CHECK_EQUAL(1,  integrals.first[idx][0]);
    BOOST_CHECK_EQUAL(10, integrals.first[idx][1]);
    BOOST_CHECK_EQUAL(1,  integrals.first[idx][2]);
    BOOST_CHECK_EQUAL(10, integrals.first[idx][3]);
    BOOST_CHECK_EQUAL(2,  integrals.first[idx][4]);
    BOOST_CHECK_EQUAL(10, integrals.first[idx][5]);
    BOOST_CHECK_EQUAL(2,  integrals.first[idx][6]);
    BOOST_CHECK_EQUAL(10, integrals.first[idx][7]);
    BOOST_CHECK_EQUAL(-1, integrals.first[idx][8]);
    BOOST_CHECK_EQUAL(-1, integrals.first[idx][9]);
    BOOST_CHECK_EQUAL(-1, integrals.first[idx][10]);
    BOOST_CHECK_EQUAL(-1, integrals.first[idx][11]);
}

/** Test for the threshold functionality of the integral parser */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_Threshold, NModeFixture)
{
    auto integrals = Vibrational::detail::NModeIntegralParser<double>(parametersFADOneBody, lattice_factory(parametersFADOneBody));
    BOOST_CHECK_EQUAL(integrals.first.size(), 78);
    BOOST_CHECK_EQUAL(integrals.second.size(), 78);
    parametersFADOneBody.set("integral_cutoff", 1.0E-6);
    integrals = Vibrational::detail::NModeIntegralParser<double>(parametersFADOneBody, lattice_factory(parametersFADOneBody));
    BOOST_CHECK_EQUAL(integrals.first.size(), 39);
    BOOST_CHECK_EQUAL(integrals.second.size(), 39);
}

/** Tests the [create_terms] method */
BOOST_FIXTURE_TEST_CASE(Test_Model_Create_Terms, NModeFixture)
{
    auto lattice = Lattice(parametersFADOneBody);
    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
    auto sizeBefore = nModeModel.hamiltonian_terms().size();
    BOOST_CHECK_EQUAL(sizeBefore, 0);
    nModeModel.create_terms();
    auto sizeAfter = nModeModel.hamiltonian_terms().size();
    BOOST_CHECK_EQUAL(sizeAfter, 78);
}

/** Checks consistency for the physical dimensions for a 1-mode system */
BOOST_FIXTURE_TEST_CASE(Test_Model_PhysDim_OneMode, NModeFixture)
{
    auto lattice = Lattice(parametersFADOneBody);
    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
    const auto& physicalDimensions0 = nModeModel.phys_dim(0);
    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), 2);
}

/** Checks consistency for the physical dimensions for a 2-mode system */
BOOST_FIXTURE_TEST_CASE(Test_Model_PhysDim_TwoMode, NModeFixture)
{
    auto lattice = Lattice(parametersFADTwoBody);
    auto nModeModel = NMode<tmatrix<double>, 2>(lattice, parametersFADTwoBody, false);
    const auto& physicalDimensions0 = nModeModel.phys_dim(0);
    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), 2);
    const auto& physicalDimensions1 = nModeModel.phys_dim(1);
    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), 2);
}

/** Checks consistency for the overall QN for a 1-mode system */
BOOST_FIXTURE_TEST_CASE(Test_Model_TotalQN_OneMode, NModeFixture)
{
    auto lattice = Lattice(parametersFADOneBody);
    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
    auto totalQN = nModeModel.total_quantum_numbers(parametersFADOneBody);
    BOOST_CHECK_EQUAL(totalQN[0], 1);
}

/** Checks consistency for the overall QN for a 1-mode system */
BOOST_FIXTURE_TEST_CASE(Test_Model_TotalQN_TwoMode, NModeFixture)
{
    auto lattice = Lattice(parametersFADTwoBody);
    auto nModeModel = NMode<tmatrix<double>, 2>(lattice, parametersFADTwoBody, false);
    auto totalQN = nModeModel.total_quantum_numbers(parametersFADTwoBody);
    BOOST_CHECK_EQUAL(totalQN[0], 1);
    BOOST_CHECK_EQUAL(totalQN[1], 1);
    // Now we use a larger lattice
    auto nModeModelWrongDim = NMode<tmatrix<double>, 3>(lattice, parametersFADTwoBody, false);
    auto totalQNWrongDim = nModeModelWrongDim.total_quantum_numbers(parametersFADTwoBody);
    BOOST_CHECK_EQUAL(totalQNWrongDim[0], 1);
    BOOST_CHECK_EQUAL(totalQNWrongDim[1], 1);
    BOOST_CHECK_EQUAL(totalQNWrongDim[2], 0);
}

/** Simple check on tags */
BOOST_FIXTURE_TEST_CASE(Test_Model_Tag_SimpleCheck_OneMode, NModeFixture)
{
    auto lattice = Lattice(parametersFADOneBody);
    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
    auto identityTag = nModeModel.identity_matrix_tag(0);
    auto fillingTag = nModeModel.filling_matrix_tag(0);
    // The nMode Hamiltonian is bosonic, so the tag should be the same
    BOOST_CHECK(identityTag == fillingTag);
}

/** Simple check on tags for the two-mode Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Model_Tag_SimpleCheck_TwoMode, NModeFixture)
{
    auto lattice = Lattice(parametersFADTwoBody);
    auto nModeModel = NMode<tmatrix<double>, 2>(lattice, parametersFADTwoBody, false);
    auto identityTag = nModeModel.filling_matrix_tag(0);
    auto fillingTag = nModeModel.filling_matrix_tag(1);
    BOOST_CHECK(identityTag != fillingTag);
}

/** Check on symbolic operator getter */
BOOST_FIXTURE_TEST_CASE(Test_Model_Symbolic_Operator_OneMode, NModeFixture)
{
    auto lattice = Lattice(parametersFADOneBody);
    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
    BOOST_CHECK(nModeModel.filling_matrix_tag(0) == nModeModel.get_operator_tag("fill", 0));
}

/** Check on symbolic operator getter for a two-mode Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Model_Symbolic_Operator_TwoMode, NModeFixture)
{
    auto lattice = Lattice(parametersFADOneBody);
    auto nModeModel = NMode<tmatrix<double>, 2>(lattice, parametersFADOneBody, false);
    BOOST_CHECK(nModeModel.filling_matrix_tag(0) == nModeModel.get_operator_tag("fill", 0));
    BOOST_CHECK(nModeModel.identity_matrix_tag(1) == nModeModel.get_operator_tag("id", 1));
}

/** Check on operator table getter */
BOOST_FIXTURE_TEST_CASE(Test_Model_Operator_Table, NModeFixture)
{
    auto lattice = Lattice(parametersFADOneBody);
    auto nModeModel = NMode<tmatrix<double>, 1>(lattice, parametersFADOneBody, false);
    auto createTag = nModeModel.get_operator_tag("bdag", 0);
    auto destroyTag = nModeModel.get_operator_tag("b", 0);
    auto countTag = nModeModel.get_operator_tag("n", 0);
    auto countTagFromTable = nModeModel.operators_table()->get_product_tag(destroyTag, createTag);
    auto siteOperatorCount = nModeModel.operators_table()->get_op(countTag);
    auto siteOperatorCountFromTable = nModeModel.operators_table()->get_op(countTagFromTable.first);
    BOOST_CHECK_CLOSE(siteOperatorCount.norm(), 1., 1.0E-16);
    BOOST_CHECK_CLOSE(siteOperatorCountFromTable.norm(), 1., 1.0E-16);
    auto differenceOperator = siteOperatorCount - siteOperatorCountFromTable;
    BOOST_CHECK_CLOSE(differenceOperator.norm(), 0., 1.0E-16);
}

#endif // HAVE_NU1
