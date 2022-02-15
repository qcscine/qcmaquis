/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2021 Laboratory for Physical Chemistry, ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@phys.ethz.ch>
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