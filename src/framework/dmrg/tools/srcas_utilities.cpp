/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#include "srcas_utilities.h"
#include "dmrg/utils/DmrgParameters.h"

#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <string>

namespace maquis {
namespace srcas {

template <typename ScalarType>  // real or complex, nmode or canonical (watson)
void SRCAS<ScalarType>::quicksort_(
    std::string dets[], ScalarType b[], int left, int right
) {
  double pivot = std::abs(b[(left + right) / 2]);
  int l = left;
  int r = right;
  while (l <= r) {
    while (std::abs(b[l]) < pivot) {
      l++;
    }
    while (std::abs(b[r]) > pivot) {
      r--;
    }
    if (l <= r) {
      ScalarType tmp = b[l];
      b[l] = b[r];
      b[r] = tmp;
      std::string ctmp = dets[l];
      dets[l] = dets[r];
      dets[r] = ctmp;
      l++;
      r--;
    }
  };
  if (left < r) {
    quicksort_(dets, b, left, r);
  }
  if (l < right) {
    quicksort_(dets, b, l, right);
  }
}

template <typename ScalarType>  // real or complex
int SRCAS<ScalarType>::getARandomOccSpinOrb_(std::vector<int> det) {
  std::vector<int> indices;
  for (int i = 1; i <= det.size(); i++) {
    if (det[i - 1] == 4) {
      indices.push_back(i);   // alpha spin orb of spatial orb i
      indices.push_back(-i);  // beta spin orb of spatial orb i
    } else if (det[i - 1] == 3) {
      indices.push_back(i);  // alpha spin orb of spatial orb i
    } else if (det[i - 1] == 2) {
      indices.push_back(-i);  // beta spin orb of spatial orb i
    }
  }
  int whichIndex = int(floor(uniformRandomNumber_() * indices.size()));
  return indices[whichIndex];
}

template <typename ScalarType>  // real or complex
int SRCAS<ScalarType>::getARandomUnoccSpinOrb_(std::vector<int> det) {
  std::vector<int> indices;
  for (int i = 1; i <= det.size(); i++) {
    if (det[i - 1] == 1) {
      indices.push_back(i);   // alpha spin orb of spatial orb i
      indices.push_back(-i);  // beta spin orb of spatial orb i
    } else if (det[i - 1] == 2) {
      indices.push_back(i);  // alpha spin orb of spatial orb i
    } else if (det[i - 1] == 3) {
      indices.push_back(-i);  // beta spin orb of spatial orb i
    }
  }
  int whichIndex = int(floor(uniformRandomNumber_() * indices.size()));
  return indices[whichIndex];
}

template <typename ScalarType>  // real or complex
bool SRCAS<ScalarType>::symmetriesFulfilled_(std::vector<int> det) {
  int nUnpaired = 0;
  int nAlpha = 0;
  int nBeta = 0;
  for (int i = 0; i < det.size(); i++) {
    if (det[i] == 4) {
      nAlpha++;
      nBeta++;
    } else if (det[i] == 3) {
      nAlpha++;
      nUnpaired++;
    } else if (det[i] == 2) {
      nBeta++;
      nUnpaired++;
    }
  }
  if (parms_["symmetry"] == "su2u1" || parms_["symmetry"] == "su2u1pg") {
    return (nUnpaired >= parms_["spin"]);
  } else {
    return (
        (nAlpha == parms_["u1_total_charge1"]) &&
        (nBeta == parms_["u1_total_charge2"])
    );
  }
}

// TODO: smarter creation of determinants (e.g. separate alpha und beta
// orbitals, create spin flipped determinant the same time)
template <typename ScalarType>  // real or complex
std::vector<int> SRCAS<ScalarType>::generateNewDet_() {
  // Start from queen
  detTmp_ = detQueen_;
  if (parms_["MODEL"] == "nmode" || parms_["MODEL"] == "watson") {
    // Loop over the modes
    for (int i = 0; i < detTmp_.size(); i++) {
      // poisson distribution centered on the current modal
      boost::poisson_distribution<> poissonDist(detTmp_[i] + 0.5);
      boost::variate_generator<boost::mt19937 &, boost::poisson_distribution<>>
          poissonRandomNumber(generator_, poissonDist);
      do {
        // Only accept a fraction of the proposed updates to stay closer to
        // reference det
        if (uniformRandomNumber_() < parms_["srcas_samplingSpeed"]) {
          detTmp_[i] = poissonRandomNumber();
        }
      } while (!(detTmp_[i] < detSpace_[i]));  // Only accept valid occupations
    }
    if (detSpace_.size() != numParticles_) {
      for (int i = 1; i < numParticles_; i++) {
        maxDetStr_ += ",";
        maxDetStr_ += parms_["Nmax"].str();
      }
      std::vector<int> tmpVec(numParticles_, std::stoi(parms_["Nmax"].str()));
      detSpace_ = std::move(tmpVec);
    }

    // TODO: better default for sampling speed
  } else if (parms_["MODEL"] == "quantum_chemistry") {
    if (parms_["symmetry"] == "su2u1" || parms_["symmetry"] == "su2u1pg") {
      throw std::runtime_error("SRCAS is does not support SU2 symmetry");
      // numParticles_ = parms_["nelec"];
    } else {
      numParticles_ =
          int(parms_["u1_total_charge1"]) + int(parms_["u1_total_charge2"]);
    }
    maxDetStr_ = "4";
    for (int i = 1; i < parms_["L"]; i++) {
      maxDetStr_ += ",4";
    }
    std::vector<int> tmpVec(parms_["L"], 4);
    detSpace_ = std::move(tmpVec);

  } else {
    throw std::runtime_error(
        "The SRCAS class supports only vibrational and electronic Hamiltonians "
        "so far"
    );
  }
}
else if (parms_["MODEL"] == "quantum_chemistry") {
  do {
    detTmp_ = detQueen_;
    // Get the number of excited electrons
    int nele_excited = geometricRandomNumber_();
    for (int i = 0; i < nele_excited; i++) {
      int annihilate = this->getARandomOccSpinOrb_(detTmp_);
      int create = this->getARandomUnoccSpinOrb_(detTmp_);
      detTmp_[abs(annihilate) - 1] -= (annihilate < 0) ? 1 : 2;
      detTmp_[abs(create) - 1] += (create < 0) ? 1 : 2;
    }
  } while (!symmetriesFulfilled_(detTmp_));  // Only accept valid occupations
}
else {
  maquis::cout << "SRCAS determinant generation NYI for calculations other "
                  "than vibrational or electronic! Abort!"
               << std::endl;
  exit(1);
}
return detTmp_;
}

template <typename ScalarType>
void SRCAS<ScalarType>::setupModel_() {
  if (parms_["MODEL"] == "nmode") {
    // Get the number of modes and the maximum occupation of each one
    numParticles_ = parms_["nmode_num_modes"];
    maxDetStr_ = parms_["nmode_num_basis"].str();
    detSpace_ = parms_["nmode_num_basis"].as<std::vector<int>>();

  } else if (parms_["MODEL"] == "watson") {
    numParticles_ = parms_["L"];
    maxDetStr_ = parms_["Nmax"].str();
    detSpace_ = parms_["Nmax"].as<std::vector<int>>();
    if (detSpace_.size() != numParticles_ && detSpace_.size() != 1) {
      throw std::runtime_error(
          "The Nmax parameter must be either a single integer, or a vector of "
          "lenght L"
      );
    }
    if (detSpace_.size() != numParticles_) {
      for (int i = 1; i < numParticles_; i++) {
        maxDetStr_ += ",";
        maxDetStr_ += parms_["Nmax"].str();
      }
      std::vector<int> tmpVec(numParticles_, std::stoi(parms_["Nmax"].str()));
      detSpace_ = std::move(tmpVec);
    }

    // TODO: better default for sampling speed
  } else if (parms_["MODEL"] == "quantum_chemistry") {
    if (parms_["symmetry"] == "su2u1" || parms_["symmetry"] == "su2u1pg") {
      throw std::runtime_error("SRCAS is does not support SU2 symmetry");
      // numParticles_ = parms_["nelec"];
    } else {
      numParticles_ =
          int(parms_["u1_total_charge1"]) + int(parms_["u1_total_charge2"]);
    }
    maxDetStr_ = "4";
    for (int i = 1; i < parms_["L"]; i++) {
      maxDetStr_ += ",4";
    }
    std::vector<int> tmpVec(parms_["L"], 4);
    detSpace_ = std::move(tmpVec);

  } else {
    throw std::runtime_error(
        "The SRCAS class supports only vibrational and electronic Hamiltonians "
        "so far"
    );
  }
}

template <typename ScalarType>
void SRCAS<ScalarType>::setupInitState_() {
  // If user set a starting det, use this, otherwise use the HF/VSCF ground
  // state
  if (parms_.is_set("init_basis_state")) {
    startingDet_ = parms_["init_basis_state"].str();
    detQueen_ = parms_["init_basis_state"].as<std::vector<int>>();

  } else {
    if (parms_["MODEL"] == "nmode" || parms_["MODEL"] == "watson") {
      startingDet_ = "0";
      for (int i = 1; i < numParticles_; i++) {
        startingDet_ += ",0";
      }
      std::vector<int> tmpVec(numParticles_, 0);
      detQueen_ = std::move(tmpVec);

    } else {
      int numDoubleOcc = numParticles_ / 2;
      startingDet_ = "";
      detQueen_.resize(parms_["L"]);
      for (int i = 0; i < numDoubleOcc; i++) {
        startingDet_ += "4,";
        detQueen_[i] = 4;
      }
      if (numParticles_ % 2) {
        startingDet_ += "3,";
        detQueen_[numDoubleOcc] = 3;
      }
      for (int i = numDoubleOcc + (numParticles_ % 2); i < parms_["L"]; i++) {
        startingDet_ += "1,";
        detQueen_[i] = 1;
      }
      startingDet_.pop_back();
    }
  }
  if (parms_["MODEL"] == "quantum_chemistry") {
    if (!symmetriesFulfilled_(detQueen_)) {
      detQueen_ = generateNewDet_();
      startingDet_ = std::to_string(detQueen_[0]);
      for (int i = 1; i < detQueen_.size(); i++) {
        startingDet_ += ",";
        startingDet_ += std::to_string(detQueen_[i]);
      }
    }
  }
  detTmp_ = detQueen_;
}

template <typename ScalarType>  // real or complex
SRCAS<ScalarType>::
    SRCAS(DmrgParameters &parameters, std::shared_ptr<InterfaceType> interface)
    : interface_(interface),
      uniformDist_(0., 1.),
      uniformRandomNumber_(generator_, uniformDist_),
      geomDist_(1.0 - parameters["srcas_samplingSpeed"]),
      geometricRandomNumber_(generator_, geomDist_),
      parms_(parameters) {
  generator_.seed(parms_["seed"]);
  this->setupModel_();
  this->setupInitState_();
}

template <typename ScalarType>
void SRCAS<ScalarType>::printSRCASSettings() {
  maquis::cout << std::endl << "----- SRCAS SETTINGS -----" << std::endl;
  maquis::cout << "MPS taken from:                             "
               << parms_["chkpfile"].str() << std::endl;
  maquis::cout << "Determinant space is:                       " << maxDetStr_
               << std::endl;
  maquis::cout << "Starting determinant is:                    " << startingDet_
               << std::endl;
  maquis::cout << "CI coeff (overlap) threshold is:            "
               << parms_["srcas_overlapThreshold"] << std::endl;
  maquis::cout << "SRCAS target completeness is:               "
               << parms_["srcas_targetCompleteness"] << std::endl;
  maquis::cout << "Maximum number of iterations is:            "
               << parms_["srcas_maxNumIterations"] << std::endl;
  maquis::cout << "Number of samples per iteration is:         "
               << parms_["srcas_numSamples"] << std::endl;
  maquis::cout << "Random number seed is:                      "
               << parms_["seed"] << std::endl;
  maquis::cout << "Sampling speed for simultaneous updates is: "
               << parms_["srcas_samplingSpeed"] << std::endl;
}

template <typename ScalarType>  // real or complex
void SRCAS<ScalarType>::run() {
  maquis::cout << std::endl
               << "----- Starting SRCAS -----" << std::endl
               << std::endl;

  // Starting det should always be added to the list
  ScalarType overlap = interface_->getCICoefficient(startingDet_);
  hashTable_[detQueen_] = overlap;
  completeness_ += pow(std::abs(overlap), 2);
  // Initialize variables that are used during the sampling
  // double sum_ci2 = 0.0;
  ScalarType ci0 = overlap;
  int nMacroIter = 0;
  int nSampled = 0;
  int nAcceptedQueen = 0;

  maquis::cout << std::setw(7) << "Iter." << std::setw(12) << "Mic. Iter."
               << std::setw(12) << "# Sample"
               << std::setw(numParticles_ * 2 + 2) << "Det." << std::setw(20)
               << "Coeff." << std::setw(20) << "Completeness" << std::endl;
  maquis::cout << std::setw(7) << nMacroIter << std::setw(12) << 0
               << std::setw(12) << nSampled << std::setw(numParticles_ * 2 + 2)
               << startingDet_ << std::setw(20) << std::setprecision(14)
               << std::fixed << overlap << std::setw(20)
               << std::setprecision(14) << std::fixed << completeness_
               << std::endl;

  do {
    // For every macroiteration generate N determinants
    for (int isample = 0; isample < parms_["srcas_numSamples"]; isample++) {
      detTmp_ = generateNewDet_();

      // Updates the data if the determinant not in map
      iter_ = hashTable_.find(detTmp_);
      if (iter_ == hashTable_.end()) {
        detTmpStr_ = std::to_string(detTmp_[0]);
        for (int i = 1; i < detTmp_.size(); i++) {
          detTmpStr_ += ",";
          detTmpStr_ += std::to_string(detTmp_[i]);
        }
        overlap = interface_->getCICoefficient(detTmpStr_);
        if (std::abs(overlap) >= parms_["srcas_overlapThreshold"]) {
          hashTable_[detTmp_] = overlap;
          nSampled++;
          completeness_ += pow(std::abs(overlap), 2.0);

          maquis::cout << std::setw(7) << nMacroIter << std::setw(12) << isample
                       << std::setw(12) << nSampled
                       << std::setw(numParticles_ * 2 + 2) << detTmpStr_
                       << std::setw(20) << std::setprecision(14) << std::fixed
                       << overlap << std::setw(20) << std::setprecision(14)
                       << std::fixed << completeness_ << std::endl;
        }
      } else {
        overlap = iter_->second;
      }

      // Determinant update (regardless of being in the hash table or not to
      // avoid getting stuck in the Markov chain) Selection criterion based on
      // CI coeff^2 in analogy to the completeness measure
      double ci_ratio = pow(std::abs(overlap), 2.0) / pow(std::abs(ci0), 2);
      double x = uniformRandomNumber_();
      if (ci_ratio > x) {
        detQueen_ = detTmp_;
        ci0 = overlap;
        nAcceptedQueen++;
      }
    }
    // sum_ci2 = calculateCompleteness_();
    nMacroIter++;

    // Prints results
    maquis::cout
        << "----------------------------------------------------------------"
        << std::endl;
    maquis::cout << "Macroiteration number:                       "
                 << nMacroIter << std::endl;
    maquis::cout << "Determinants sampled above the CI threshold: " << nSampled
                 << std::endl;
    maquis::cout << "Determinants accepted as queens:             "
                 << nAcceptedQueen << std::endl;
    maquis::cout << "Current completeness (\\sum(ci^2)):           "
                 << completeness_ << std::endl;

  } while ((completeness_ < parms_["srcas_targetCompleteness"]) &&
           (nMacroIter < parms_["srcas_maxNumIterations"]));
  // Final completeness
  // completeness_ = sum_ci2;
}

// +---------------+
//   FINAL PRINTING
// +---------------+
template <typename ScalarType>  // real or complex, nmode or canonical (watson)
void SRCAS<ScalarType>::printResults() {
  maquis::cout
      << "----------------------------------------------------------------"
      << std::endl;
  maquis::cout << std::endl << "--- Finished SRCAS ---" << std::endl;
  maquis::cout << "Final completeness is:                " << completeness_
               << std::endl;
  maquis::cout << "# of stored determinants is:          " << hashTable_.size()
               << std::endl;

  ScalarType CIs_show[hashTable_.size()];    // CI value
  std::string dets_show[hashTable_.size()];  // dets represent
  int i = 0;
  int det_length = hashTable_.begin()->first.size();
  maquis::cout
      << std::endl
      << "-------------DETERMINANTS ABOVE OVERLAP THRESHOLD--------------------"
      << std::endl
      << std::endl;
  for (iter_ = hashTable_.begin(); iter_ != hashTable_.end(); iter_++) {
    // Local initialization that is later used for sorting
    std::string ctmp;
    CIs_show[i] = iter_->second;
    for (int p = 0; p < det_length; p++) {
      ctmp = ctmp + boost::lexical_cast<std::string>(iter_->first[p]);
    }
    dets_show[i] = ctmp;
    i++;
  }
  // Final sorting
  quicksort_(dets_show, CIs_show, 0, hashTable_.size() - 1);
  // Output the entire ordered list
  maquis::cout << std::fixed << std::setprecision(10);
  for (int i = 0; i < hashTable_.size(); i++) {
    maquis::cout << " Determinant " << dets_show[hashTable_.size() - i - 1]
                 << " with ";
    // if (CIs_show[hashTable_.size()-i-1]>0) maquis::cout << " ";
    maquis::cout << CIs_show[hashTable_.size() - i - 1] << " is number "
                 << i + 1 << std::endl;
  }
}

template <typename ScalarType>  // real or complex, nmode or canonical (watson)
std::vector<int> SRCAS<ScalarType>::getCurrentQueen() {
  return detQueen_;
}

template <typename ScalarType>  // real or complex, nmode or canonical (watson)
std::map<std::vector<int>, ScalarType> SRCAS<ScalarType>::getDetTable() {
  return hashTable_;
}

template <typename ScalarType>  // real or complex, nmode or canonical (watson)
double SRCAS<ScalarType>::getCompleteness() {
  return completeness_;
}

// Explicit template instantiation
template class SRCAS<double>;
template class SRCAS<std::complex<double>>;
}  // namespace srcas
}  // namespace maquis
