/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef VIBRATIONAL_MODEL_TRAIT_CLASS_H
#define VIBRATIONAL_MODEL_TRAIT_CLASS_H

#include "dmrg/block_matrix/symmetry.h"

/** @brief Trait class containing data required for vibrational models */
template<class SymmGroup>
class VibrationalModelTraitClass {};

/**
 * @brief Class specialization for the TrivialGroup
 * The trait class includes the maximum number of coupling terms that are supported.
 */
template<>
class VibrationalModelTraitClass<TrivialGroup> {
public:
  static constexpr int maximumNumberOfCouplings = 6;
};

/**
 * @brief Class specialization for the TrivialGroup
 * The trait class includes the maximum number of coupling terms that are supported.
 */
template<>
class VibrationalModelTraitClass<U1> {
public:
  static constexpr int maximumNumberOfCouplings = 2;
};

#endif // VIBRATIONAL_MODEL_TRAIT_CLASS_H
