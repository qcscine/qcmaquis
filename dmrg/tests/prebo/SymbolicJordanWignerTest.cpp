/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <alberto.baiardi@phys.chem.ethz.ch>
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

#define BOOST_TEST_MODULE SymbolicJordanWigner

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/coded/nu1_SymbolicJordanWigner.hpp"

/** Checks that the constructor works */
BOOST_AUTO_TEST_CASE( PreBO_SymbolicJW_BosonicOperator )
{
    auto createOp = SymbolicOperator(0, OpType::Create, 0, Spin::None);
    auto destroyOp = SymbolicOperator(1, OpType::Annihilate, 0, Spin::None);
    auto vectorOfOperators = std::vector<SymbolicOperator>{createOp, destroyOp};
    auto jwSingleExcitationBoson = SymbolicJordanWigner(vectorOfOperators);
    BOOST_CHECK_EQUAL(vectorOfOperators.size(), 2);
}

/** 
 * Checks the JW transformation for a one-body excitation for fermionic systems 
 * Note that for alpha excitation the filling operator is attached to the 
 * operator with the *lowest* index.
 */
BOOST_AUTO_TEST_CASE( PreBO_SymbolicJW_FermionicOneBodyExcitation )
{
    auto createOpF = SymbolicOperator(0, OpType::Create, 0, Spin::Up);
    auto destroyOpF = SymbolicOperator(2, OpType::Annihilate, 0, Spin::Up);
    auto vectorOfOperators = std::vector<SymbolicOperator>{createOpF, destroyOpF};
    auto jwSingleExcitationFermion = SymbolicJordanWigner(vectorOfOperators);
    // Preliminary check
    BOOST_CHECK_EQUAL(jwSingleExcitationFermion.getSymOpStr().size(), 3);
    // Checks that the second operator is a filling operator
    bool checkFilling = jwSingleExcitationFermion.getSymOpStr()[1].getOpType() == OpType::Filling;
    BOOST_CHECK(checkFilling);
}

/** 
 * Checks the JW transformation for the fully-diagonal two-body excitation
 */
BOOST_AUTO_TEST_CASE( PreBO_SymbolicJW_FermionicTwoBodyExcitationSameSite )
{
    auto createOpF1 = SymbolicOperator(0, OpType::Create, 0, Spin::Up);
    auto createOpF2 = SymbolicOperator(0, OpType::Create, 0, Spin::Down);
    auto destroyOpF2 = SymbolicOperator(0, OpType::Annihilate, 0, Spin::Down);
    auto destroyOpF1 = SymbolicOperator(0, OpType::Annihilate, 0, Spin::Up);
    auto vectorOfOperators = std::vector<SymbolicOperator>{createOpF1, createOpF2, destroyOpF2, destroyOpF1};
    auto jwSingleExcitationFermion = SymbolicJordanWigner(vectorOfOperators);
    auto operatorVector = jwSingleExcitationFermion.getSymOpStr();
    // Preliminary check
    BOOST_CHECK_EQUAL(operatorVector.size(), 4);
    // Checks that the second operator is a filling operator
    for (int iOp = 0; iOp < 4; iOp++) {
        bool checkOpConsistent = operatorVector[iOp].getOpType() == vectorOfOperators[iOp].getOpType();
        BOOST_CHECK(checkOpConsistent);
    }
}

/** 
 * Checks the JW transformation for a semi-diagonal two-body operator.
 */
BOOST_AUTO_TEST_CASE( PreBO_SymbolicJW_FermionicTwoBodyExcitationSemidiagonal )
{
    auto createFermion = SymbolicOperator(0, OpType::Create, 1, Spin::Up);
    auto destroyFermion = SymbolicOperator(1, OpType::Annihilate, 1, Spin::Up);
    auto createBoson = SymbolicOperator(2, OpType::Create, 1, Spin::Down);
    auto destroyBoson = SymbolicOperator(3, OpType::Annihilate, 1, Spin::Down);
    auto vectorOfOperators = std::vector<SymbolicOperator>{createFermion, destroyFermion, createBoson, destroyBoson};
    auto jwSingleExcitationFermion = SymbolicJordanWigner(vectorOfOperators);
    auto operatorVector = jwSingleExcitationFermion.getSymOpStr();
    BOOST_CHECK_EQUAL(operatorVector.size(), 6);
    bool checkFilling = jwSingleExcitationFermion.getSymOpStr()[1].getOpType() == OpType::Filling;
    BOOST_CHECK(checkFilling);
    bool checkPos = jwSingleExcitationFermion.getSymOpStr()[1].getSite() == 0;
    BOOST_CHECK(checkPos);
}