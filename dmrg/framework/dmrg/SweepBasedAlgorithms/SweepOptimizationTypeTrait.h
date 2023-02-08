/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <alberto.baiardi@phys.chem.ethz.ch>
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

#ifndef SWEEP_OPTIMIZATION_TYPE_TRAIT
#define SWEEP_OPTIMIZATION_TYPE_TRAIT

#include "dmrg/mp_tensors/mpstensor.h"

enum class SweepOptimizationType {SingleSite, TwoSite};

enum class SweepDirectionType {Forward, Backward, EndOfLattice};

enum class GrowBoundaryModality {LeftToRight, RightToLeft};

/** @brief Declaration of the trait class */
template<SweepOptimizationType SweepType>
class SweepOptimizationTypeTrait {};

/** @brief Specialization of the class for the SS case */
template<>
class SweepOptimizationTypeTrait<SweepOptimizationType::SingleSite> {
public:
  /** @brief Gets the index of the left boundary for a given site */
  static int getIndexOfLeftBoundary(int L, int microiteration) {
    return convertMicroIterationToSite(L, microiteration);
  }

  /** @brief Gets the index of the right boundary for a given site */
  static int getIndexOfRightBoundary(int L, int microiteration) {
    return convertMicroIterationToSite(L, microiteration)+1;
  }

  /**
   * @brief Gets the number of microiterations for a given
   * Note that the case in which L == 1/2 requires a special treatment.
   */
  static int getNumberOfMicroiterations(int L) {
    int numberOfMicroiterations;
    if (L == 1)
      numberOfMicroiterations = 1;
    else if (L == 2)
      numberOfMicroiterations = 2;
    else
      numberOfMicroiterations = 2*(L-1);
    return numberOfMicroiterations;
  }

  /** @brief Gets the direction of a sweep */
  static SweepDirectionType getSweepDirection(int L, int i) {
    SweepDirectionType ret;
    if (i < L)
      ret = SweepDirectionType::Forward;
    else if (i < getNumberOfMicroiterations(L))
      ret = SweepDirectionType::Backward;
    else
      ret = SweepDirectionType::EndOfLattice;
    return ret;
    // return (i < L) ? SweepDirectionType::Forward : SweepDirectionType::Backward;
  }

  /** @brief Simple function converting microiteration index to the optimization site */
  static int convertMicroIterationToSite(int L, int i) {
    assert(i < getNumberOfMicroiterations(L));
    return (i < L) ? i : 2*L-2-i;
  }

  /** @brief Returns true if at the next microiteration the direction will be reversed */
  static bool changeDirectionNextMicroiteration(int L, int i) {
    return (i == L-1) || (i == getNumberOfMicroiterations(L)-1);
  }

  /** @brief Gets the string identifier of the simulation type */
  static std::string getSimulationTypeName() {
    return std::string("Single site");
  }

  // Static members
  static const bool countEndSiteTwice_=false;
};

/** @brief Specialization of the class for the TS case */
template<>
class SweepOptimizationTypeTrait<SweepOptimizationType::TwoSite> {
public:
  /** @brief Gets the index of the left boundary for a given site */
  static inline int getIndexOfLeftBoundary(int L, int microiteration) {
    return convertMicroIterationToSite(L, microiteration);
  }

  /** @brief Gets the index of the right boundary for a given site */
  static inline int getIndexOfRightBoundary(int L, int microiteration) {
    return convertMicroIterationToSite(L, microiteration)+2;
  }

  /** @brief Simple function converting microiteration index to the optimization site */
  static int convertMicroIterationToSite(int L, int i) {
    assert(i < getNumberOfMicroiterations(L));
    return (i < L-1) ? i : 2*L-4 - i;
  }

  /** @brief Gets the number of microiterations for a given */
  static int getNumberOfMicroiterations(int L) {
    int numberOfMicroiterations;
    if (L == 2)
      numberOfMicroiterations = 1;
    else if (L == 1)
      throw std::runtime_error("TS optimizer not available for 1-site lattices");
    else
      numberOfMicroiterations = 2*(L-2);
    return numberOfMicroiterations;
  }

  /** @brief Gets the direction of a sweep */
  static SweepDirectionType getSweepDirection(int L, int i) {
    SweepDirectionType ret;
    if (i < L-1)
      ret = SweepDirectionType::Forward;
    else if (i < getNumberOfMicroiterations(L))
      ret = SweepDirectionType::Backward;
    else
      ret = SweepDirectionType::EndOfLattice;
    return ret;
    // return (i < L-1) ? SweepDirectionType::Forward : SweepDirectionType::Backward;
  }

  /** @brief Returns true if at the next microiteration the direction will be reversed */
  static bool changeDirectionNextMicroiteration(int L, int i) {
    return (i == L-2) || (i == getNumberOfMicroiterations(L)-1);
  }

  /** @brief Gets the string identifier of the simulation type */
  static std::string getSimulationTypeName() {
    return std::string("Two site");
  }

  // Static members
  static const bool countEndSiteTwice_=true;
};

#endif // SWEEP_OPTIMIZATION_TYPE_TRAIT
