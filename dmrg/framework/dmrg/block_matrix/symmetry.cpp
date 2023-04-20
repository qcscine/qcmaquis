/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "utils/io.hpp" // has to be first include because of impi
#include "dmrg/block_matrix/symmetry.h"

const TrivialGroup::charge  TrivialGroup::IdentityCharge;
const U1::charge            U1::IdentityCharge;
const Ztwo::charge          Ztwo::IdentityCharge;
