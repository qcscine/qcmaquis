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