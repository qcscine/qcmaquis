/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE SymbolicOperator

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/prebo/nu1/model.hpp"
#include "dmrg/models/prebo/nu1/nu1_SymbolicJordanWigner.hpp"
#include <iostream>

/** Checks that the constructor works */
BOOST_AUTO_TEST_CASE( PreBO_SymbolicOperator_CheckConstructor )
{
    auto createOp = SymbolicOperator(0, OpType::Create, 0);
    BOOST_CHECK_EQUAL(createOp.getSite(), 0);
}

/** Checks that the copy constructor works */
BOOST_AUTO_TEST_CASE( PreBO_SymbolicOperator_CheckCopyConstructor )
{
    auto createOp = SymbolicOperator(0, OpType::Create, 1);
    auto copiedOp = SymbolicOperator(createOp);
    BOOST_CHECK_EQUAL(copiedOp.getPartType(), 1);
}

/** Checks that the modified copy constructor works */
BOOST_AUTO_TEST_CASE( PreBO_SymbolicOperator_CheckModifiedCopyConstructor )
{
    auto createOp = SymbolicOperator(0, OpType::Create, 1);
    auto copiedOp = SymbolicOperator(createOp, OpType::Ident);
    bool valid = copiedOp.getSpin() == Spin::None;
    BOOST_CHECK(valid);
}