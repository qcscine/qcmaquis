/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2023 Reiher Group, ETH Zurich
 *               2023- by Nina Glaser <nglaser@phys.chem.ethz.ch>
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

#ifndef SRCAS_UTILITIES_H
#define SRCAS_UTILITIES_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

#include <boost/random.hpp>
#include <memory>
#include <string>

namespace maquis {
namespace srcas {

/**
 * @brief SRCAS class responsible for the sampling of an MPS
 *
 * @tparam ScalarType double or complex
 */
template <typename ScalarType> // real or complex
class SRCAS {
  using InterfaceType = maquis::DMRGInterface<ScalarType>;

public:
  /**
   * @brief Constructor
   *
   * @param parameters qcmaquis parameter object
   * @param interface qcmaquis dmrg interface
   */
  SRCAS(DmrgParameters &parameters, std::shared_ptr<InterfaceType> interface);
  /** @brief Run SRCAS sampling */
  void run();
  /** @brief Print SRCAS settings */
  void printSRCASSettings();
  /** @brief Print SRCAS results */
  void printResults();
  /** @brief get last determinant queen */
  std::vector<int> getCurrentQueen();
  /** @brief get map with all sampled determinants above the threshold */
  std::map<std::vector<int>, ScalarType> getDetTable();
  /**
   * @brief get achieved completeness
   *
   * @return double the completeness
   */
  double getCompleteness();

private:
  /** @brief generate a new determinant form current queen */
  std::vector<int> generateNewDet();

  /** @brief Evaluate the current completeness */
  double calculateCompleteness();
  void quicksort(std::string dets[], ScalarType b[], int left, int right);

  // For electronic case
  int getARandomOccSpinOrb(std::vector<int> det);
  int getARandomUnoccSpinOrb(std::vector<int> det);
  bool symmetriesFulfilled(std::vector<int> det);

  boost::mt19937 generator_;
  boost::uniform_real<> uniformDist_;
  boost::geometric_distribution<double> geomDist_;
  boost::variate_generator<boost::mt19937 &, boost::uniform_real<double>>
      uniformRandomNumber_;
  boost::variate_generator<boost::mt19937 &,
                           boost::geometric_distribution<double>>
      geometricRandomNumber_;

  DmrgParameters &parms_;
  std::shared_ptr<InterfaceType> interface_;

  std::string startingDet_, maxDetStr_, detTmpStr_;
  std::vector<int> detQueen_, detTmp_, detSpace_;
  int numParticles_; // this is either modes or electrons, for the vibrational
                     // or the electronic case, respectively
  double completeness_;

  std::map<std::vector<int>, ScalarType> hashTable_;
  typename std::map<std::vector<int>, ScalarType>::iterator iter_;
};
} // namespace srcas
} // namespace maquis

#endif
