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

enum class SweepDirectionType {Forward, Backward};

/** @brief Declaration of the trait class */
template<SweepOptimizationType SweepType>
class SweepOptimizationTypeTrait {};

/** @brief Specialization of the class for the SS case */
template<>
class SweepOptimizationTypeTrait<SweepOptimizationType::SingleSite> {
public:
  /** @brief Gets the index of the left boundary for a given site */
  static int getIndexOfLeftBoundary(int site, SweepDirectionType sweepDirection) { return site; }
  
  /** @brief Gets the index of the right boundary for a given site */
  static int getIndexOfRightBoundary(int site, SweepDirectionType sweepDirection) { return site+1; }

  /** @brief Gets the index of the left boundary associated with the following sweep for a given site */
  static int getIndexOfNextLeftBoundary(int site, SweepDirectionType sweepDirection) {
    return (sweepDirection == SweepDirectionType::Forward) ? getIndexOfLeftBoundary(site+1, sweepDirection) 
                                                           : getIndexOfLeftBoundary(site-1, sweepDirection);
  }
  
  /** @brief Gets the index of the right boundary associated with the following sweep for a given site */
  static int getIndexOfNextRightBoundary(int site, SweepDirectionType sweepDirection) {
    return (sweepDirection == SweepDirectionType::Forward) ? getIndexOfLeftBoundary(site+1, sweepDirection) 
                                                           : getIndexOfLeftBoundary(site-1, sweepDirection);
  }

  /** @brief Gets the upper boundary of the loop over the microiterations */
  static int getLastSite(int L) {
    return L;
  }

  // Static members
  static const bool countEndSiteTwice_=false;
};

/** @brief Specialization of the class for the TS case */
template<>
class SweepOptimizationTypeTrait<SweepOptimizationType::TwoSite> {
public:
  /** @brief Gets the index of the left boundary for a given site */
  static inline int getIndexOfLeftBoundary(int site, SweepDirectionType sweepDirection) { 
    return (sweepDirection == SweepDirectionType::Forward) ? site : site-1; 
  }
  
  /** @brief Gets the index of the right boundary for a given site */
  static inline int getIndexOfRightBoundary(int site, SweepDirectionType sweepDirection) { 
    return (sweepDirection == SweepDirectionType::Forward) ? site+2 : site+1;
  }

  /** @brief Gets the index of the left boundary associated with the following sweep for a given site */
  static int getIndexOfNextLeftBoundary(int site, SweepDirectionType sweepDirection) {
    return (sweepDirection == SweepDirectionType::Forward) ? getIndexOfLeftBoundary(site+1, sweepDirection) 
                                                           : getIndexOfLeftBoundary(site-1, sweepDirection);
  }
  
  /** @brief Gets the index of the right boundary associated with the following sweep for a given site */
  static int getIndexOfNextRightBoundary(int site, SweepDirectionType sweepDirection) {
    return (sweepDirection == SweepDirectionType::Forward) ? getIndexOfLeftBoundary(site+1, sweepDirection) 
                                                           : getIndexOfLeftBoundary(site-1, sweepDirection);
  }

  /** @brief Gets the upper boundary of the loop over the microiterations */
  static int getLastSite(int L) {
    return L-1;
  }

  // Static members
  static const bool countEndSiteTwice_=true;
};

#endif // SWEEP_OPTIMIZATION_TYPE_TRAIT