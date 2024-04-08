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


// TODO: change key of map for python bindings
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
  /**
   * @brief get last determinant queen 
   *
   * @return vector representation of current queen 
   **/ 
  std::vector<int> getCurrentQueen();
  /**
   * @brief getter map with all sampled determinants above the threshold 
   *
   * @retrun map with vector representation as key and correpsonding value
   **/
  std::map<std::vector<int>, ScalarType> getDetTable();
  /**
   * @brief get achieved completeness
   *
   * @return double the completeness
   */
  double getCompleteness();

private:
  /** @brief setup SRCAS for different models */
  void setupModel();
  /** @brief generate a new determinant form current queen */
  std::vector<int> generateNewDet_();
  /** @brief Evaluate the current completeness */
  double calculateCompleteness_();
  /** @brief a simple quicksort 
   *
   * @param dets 
   * @param b 
   * @param left  
   * @param right
   **/
  void quicksort_(std::string dets[], ScalarType b[], int left, int right);
  /** @brief get a random occupied orbital for det
   *
   * @param det vector representation of a Determinant
   **/
  int getARandomOccSpinOrb_(std::vector<int> det);
  /** @brief get a random virtual orbital for det
   *
   * @param det vector representation of a Determinant
   **/
  int getARandomUnoccSpinOrb_(std::vector<int> det);
  /** @brief get a random virtual orbital for det
   *
   * @param det vector representation of a Determinant
   **/
  bool symmetriesFulfilled_(std::vector<int> det);
  /** @brief boost random number generator **/
  boost::mt19937 generator_;
  /** @brief boost uniform distribution **/
  boost::uniform_real<> uniformDist_;
  /** @brief boost geometric distribution **/
  boost::geometric_distribution<double> geomDist_;
  /** @brief boost uniform distribution generator **/
  boost::variate_generator<boost::mt19937 &, boost::uniform_real<double>> uniformRandomNumber_;
  /** @brief boost geometric distribution generator **/
  boost::variate_generator<boost::mt19937 &, boost::geometric_distribution<double>> geometricRandomNumber_;
  /** @brief all DMRG parameters **/
  DmrgParameters &parms_;
  /** @brief DMRG interface **/
  std::shared_ptr<InterfaceType> interface_;
  /** @brief string representation of initial determinant **/
  std::string startingDet_;
  /** @brief string representation of maximal determinant ?? **/
  std::string maxDetStr_;
  /** @brief string representation of current determinant **/
  std::string detTmpStr_;
  /** @brief current queen **/
  std::vector<int> detQueen_;
  /** @brief current determinant **/
  std::vector<int> detTmp_;
  /** @brief determines the possible space for each determinant **/
  std::vector<int> detSpace_;
  /** @brief modes or electrons **/
  int numParticles_; 
  /** @brief current completeness **/
  double completeness_;
  /** @brief all sampled coeffs above threshold **/
  std::map<std::vector<int>, ScalarType> hashTable_;
  
  typename std::map<std::vector<int>, ScalarType>::iterator iter_;

};
} // namespace srcas
} // namespace maquis

#endif
