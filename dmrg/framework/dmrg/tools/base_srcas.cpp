#include "base_srcas.h"
#include "determinant.h"
#include "vib_onv.h"
#include <cmath>

namespace maquis {
namespace srcas {

template<typename ScalarType, class T>
BaseSRCAS<ScalarType, T>::BaseSRCAS(DmrgParameters& parameters, std::shared_ptr<InterfaceType> interface)
  : interface_(interface),
    uniformDist_(0., 1.),
    uniformRandomNumber_(generator_, uniformDist_),
    geomDist_(1.0 - parameters["srcas_samplingSpeed"]),
    geometricRandomNumber_(generator_, geomDist_),
    parms_(parameters) {
  generator_.seed(parms_["seed"]);
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::evaluateSpecificONVs_(const std::vector<std::string>& extra_onvs) {
  for (const auto& onv : extra_onvs) {
    // ensure to not double count onvs
    iter_ = hashTable_.find(onv);
    if (iter_ != hashTable_.end()) {
      continue;
    }

    ScalarType overlap = interface_->getCICoefficient(onv);
    hashTable_[onv] = overlap;
    completeness_ += std::pow(std::abs(overlap), 2);

    maquis::cout << std::setw(12) << "extra" << std::setw(12) << "" << std::setw(parms_["L"] * 2 + 2) << onv
                 << std::setw(20) << std::setprecision(14) << std::fixed << overlap << std::setw(20)
                 << std::setprecision(14) << std::fixed << completeness_ << std::endl;
  }
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::run(const std::vector<std::string>& extra_onvs) {
  this->evaluateInitial_();
  this->evaluateSpecificONVs_(extra_onvs);
  this->sample_();
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::evaluateInitial_() {
  maquis::cout << "\n----- Starting SRCAS -----\n\n";
  // Starting onv should always be added to the list
  auto queenString = queen_.string();

  maquis::cout << "    queen: " << queenString << std::endl;

  ScalarType overlap = interface_->getCICoefficient(queenString);
  hashTable_[queenString] = overlap;
  completeness_ += std::pow(std::abs(overlap), 2);

  maquis::cout << std::setw(12) << "Iteration" << std::setw(12) << "# Sampled" << std::setw(parms_["L"] * 2 + 2)
               << "ONV " << std::setw(20) << "CI Coefficient" << std::setw(20) << "Completeness" << std::endl;
  maquis::cout << std::setw(12) << "initial" << std::setw(12) << "" << std::setw(parms_["L"] * 2 + 2) << queen_.string()
               << std::setw(20) << std::setprecision(14) << std::fixed << overlap << std::setw(20)
               << std::setprecision(14) << std::fixed << completeness_ << std::endl;
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::sample_() {
  ScalarType ci0 = 0;
  int nSampled = 0;

  // For every macroiteration generate N occupation number vectors
  for (int isample = 1; isample <= parms_["srcas_numSamples"]; isample++) {
    ScalarType overlap = 0;
    auto tmpONV = generateNewONV_();

    // Updates the data if the onv not in map
    iter_ = hashTable_.find(tmpONV.string());
    if (iter_ == hashTable_.end()) {
      overlap = interface_->getCICoefficient(tmpONV.string());
      if (std::abs(overlap) >= parms_["srcas_overlapThreshold"]) {
        hashTable_[tmpONV.string()] = overlap;
        nSampled++;
        completeness_ += std::pow(std::abs(overlap), 2.0);

        maquis::cout << std::setw(12) << isample << std::setw(12) << nSampled << std::setw(parms_["L"] * 2 + 2)
                     << tmpONV.string() << std::setw(20) << std::setprecision(14) << std::fixed << overlap
                     << std::setw(20) << std::setprecision(14) << std::fixed << completeness_ << std::endl;

        try {
          auto otherONV = generateSymmetricDeterminant_(tmpONV);
          if (!(otherONV == tmpONV)) {
            ScalarType newOverlap = interface_->getCICoefficient(otherONV.string());
            if (std::abs(newOverlap) >= parms_["srcas_overlapThreshold"]) {
              hashTable_[otherONV.string()] = newOverlap;
              completeness_ += std::pow(std::abs(newOverlap), 2.0);

              maquis::cout << std::setw(12) << "spinflip" << std::setw(12) << nSampled << std::setw(parms_["L"] * 2 + 2)
                           << otherONV.string() << std::setw(20) << std::setprecision(14) << std::fixed << newOverlap
                           << std::setw(20) << std::setprecision(14) << std::fixed << completeness_ << std::endl;
            }
          }
        }
        catch (...) {
        }
      }
    }
    else {
      overlap = iter_->second;
    }

    // ONV update (regardless of being in the hash table or not to
    // avoid getting stuck in the Markov chain) Selection criterion based on
    // CI coeff^2 in analogy to the completeness measure
    double ci_ratio = pow(std::abs(overlap), 2.0) / pow(std::abs(ci0), 2);
    if (ci_ratio > uniformRandomNumber_()) {
      queen_ = tmpONV;
      ci0 = overlap;
      maquis::cout << "new queen: " << queen_.string() << std::endl;
    }
    // }

    if (completeness_ > parms_["srcas_targetCompleteness"]) {
      maquis::cout << "SRCAS reached target completeness of " << parms_["srcas_targetCompleteness"] << std::endl;
      maquis::cout << "SRCAS current completeness           " << completeness_ << std::endl;
      break;
    }
  }
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::run() {
  this->evaluateInitial_();
  this->sample_();
}
/**
 * @brief get last ONV queen
 *
 * @return vector representation of current queen
 **/
template<typename ScalarType, class T>
std::vector<int> BaseSRCAS<ScalarType, T>::currentQueen() const {
  return queen_.vector();
}
/**
 * @brief getter map with all sampled ONV above the threshold
 *
 * @retrun map with vector representation as key and correpsonding value
 **/
template<typename ScalarType, class T>
const std::map<std::string, ScalarType>& BaseSRCAS<ScalarType, T>::sampledTable() const {
  return hashTable_;
}
/**
 * @brief get achieved completeness
 *
 * @return double the completeness
 */
template<typename ScalarType, class T>
double BaseSRCAS<ScalarType, T>::completeness() const {
  return completeness_;
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::printSettings() const {
  maquis::cout << "\n----- SRCAS settings -----\n";
  maquis::cout << "MPS checkpoint:                          " << parms_["chkpfile"].str() << "\n";
  maquis::cout << "Starting ONV:                            " << queen_.string() << "\n";
  maquis::cout << "CI coeff (overlap) threshold:            " << parms_["srcas_overlapThreshold"] << "\n";
  maquis::cout << "SRCAS target completeness:               " << parms_["srcas_targetCompleteness"] << "\n";
  maquis::cout << "Maximum number of iterations:            " << parms_["srcas_maxNumIterations"] << "\n";
  maquis::cout << "Number of samples per iteration:         " << parms_["srcas_numSamples"] << "\n";
  maquis::cout << "Random number seed:                      " << parms_["seed"] << "\n";
  maquis::cout << "Sampling speed for simultaneous updates: " << parms_["srcas_samplingSpeed"] << "\n";
  printSettingImpl_();
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::printSettingImpl_() const {
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::printResults() const {
  maquis::cout << "\n----- SRCAS summary -----\n";
  maquis::cout << "Final completeness: " << completeness_ << "\n";
  maquis::cout << "# of stored ONVs:   " << hashTable_.size() << "\n";

  maquis::cout << "\n----- ONVs above overlap threshold -----\n\n";

  std::vector<std::pair<std::string, ScalarType>> onvValuePairs;
  for (auto itr = hashTable_.begin(); itr != hashTable_.end(); itr++) {
    onvValuePairs.push_back(*itr);
  }
  std::sort(onvValuePairs.begin(), onvValuePairs.end(),
            [=](std::pair<std::string, ScalarType>& a, std::pair<std::string, ScalarType>& b) {
              return std::abs(a.second) > std::abs(b.second);
            });

  maquis::cout << std::setw(10) << "index" << std::setw(parms_["L"] * 2 + 2) << "ONV " << std::setw(20)
               << "CI Coefficient" << std::setw(20) << "Completeness" << std::endl;
  ScalarType tmpCompleteness = 0.0;
  int count = 0;
  for (const auto& i : onvValuePairs) {
    count++;
    tmpCompleteness += std::pow(i.second, 2);
    maquis::cout << std::setw(10) << count << std::setw(parms_["L"] * 2 + 2) << i.first << std::setw(20)
                 << std::setprecision(14) << std::fixed << i.second << std::setw(20) << std::setprecision(14)
                 << std::fixed << tmpCompleteness << std::endl;
  }
}

template class BaseSRCAS<double, Determinant>;
template class BaseSRCAS<std::complex<double>, Determinant>;
template class BaseSRCAS<double, VibONV>;
template class BaseSRCAS<std::complex<double>, VibONV>;

} // namespace srcas
} // namespace maquis
