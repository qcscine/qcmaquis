/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#ifdef USE_AMBIENT
#include "dmrg/mp_tensors/contractions/detail/ambient.hpp"
#else
//#include "dmrg/mp_tensors/contractions/impl/alps.hpp"
#include "dmrg/mp_tensors/contractions/detail/memsave.hpp"
#endif

#include "dmrg/mp_tensors/contractions/abelian/special.hpp"

#include "dmrg/mp_tensors/contractions/engine.h"

#endif
