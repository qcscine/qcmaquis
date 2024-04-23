
#include "base_srcas.h"
#include "determinant.h"
#include "vib_onv.h"

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
void BaseSRCAS<ScalarType, T>::run() {
  maquis::cout << "\n----- Starting SRCAS -----\n\n";
  // Starting onv should always be added to the list
  auto queenString = queen_.string();
  maquis::cout << "queen: " << queenString << std::endl;
  ScalarType overlap = interface_->getCICoefficient(queenString);
  hashTable_[queenString] = overlap;
  completeness_ += std::pow(std::abs(overlap), 2);

  ScalarType ci0 = overlap;
  int nSampled = 0;
  int nAcceptedQueen = 0;

  maquis::cout << std::setw(12) << "Mic. Iter." << std::setw(12) << "# Sample" << std::setw(parms_["L"] * 2 + 2)
               << "ONV " << std::setw(20) << "Coeff." << std::setw(20) << "Completeness"
               << "\n";
  maquis::cout << std::setw(12) << 1 << std::setw(12) << nSampled << std::setw(parms_["L"] * 2 + 2) << queen_.string()
               << std::setw(20) << std::setprecision(14) << std::fixed << overlap << std::setw(20)
               << std::setprecision(14) << std::fixed << completeness_ << std::endl;

  // For every macroiteration generate N occupation number vectors
  for (int isample = 2; isample <= parms_["srcas_numSamples"]; isample++) {
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
      nAcceptedQueen++;
      maquis::cout << "new queen: " << queen_.string() << std::endl;
    }
    // }

    if (completeness_ > parms_["srcas_targetCompleteness"]) {
      maquis::cout << "SRCAS reached target completeness of " << parms_["srcas_targetCompleteness"] << "\n";
      maquis::cout << "SRCAS current completeness           " << completeness_ << "\n";
      break;
    }
  }
}
/**
 * @brief get last ONV queen
 *
 * @return vector representation of current queen
 **/
template<typename ScalarType, class T>
T BaseSRCAS<ScalarType, T>::currentQueen() const {
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
  maquis::cout << "\n----- SRCAS SETTINGS -----\n";
  maquis::cout << "MPS taken from:                             " << parms_["chkpfile"].str() << "\n";
  maquis::cout << "Starting ONV is:                            " << queen_.string() << "\n";
  maquis::cout << "CI coeff (overlap) threshold is:            " << parms_["srcas_overlapThreshold"] << "\n";
  maquis::cout << "SRCAS target completeness is:               " << parms_["srcas_targetCompleteness"] << "\n";
  maquis::cout << "Maximum number of iterations is:            " << parms_["srcas_maxNumIterations"] << "\n";
  maquis::cout << "Number of samples per iteration is:         " << parms_["srcas_numSamples"] << "\n";
  maquis::cout << "Random number seed is:                      " << parms_["seed"] << "\n";
  maquis::cout << "Sampling speed for simultaneous updates is: " << parms_["srcas_samplingSpeed"] << "\n";
  printSettingImpl_();
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::printSettingImpl_() const {
}

template<typename ScalarType, class T>
void BaseSRCAS<ScalarType, T>::printResults() const {
  maquis::cout << "----------------------------------------------------------------\n";
  maquis::cout << "\n--- Finished SRCAS ---\n";
  maquis::cout << "Final completeness is:                " << completeness_ << "\n";
  maquis::cout << "# of stored ONVs is:                  " << hashTable_.size() << "\n";

  maquis::cout << "\n--------------ONVs ABOVE OVERLAP THRESHOLD--------------------\n\n";

  std::vector<std::pair<std::string, ScalarType>> onvValuePairs;
  for (auto itr = hashTable_.begin(); itr != hashTable_.end(); itr++) {
    onvValuePairs.push_back(*itr);
  }
  std::sort(onvValuePairs.begin(), onvValuePairs.end(),
            [=](std::pair<std::string, ScalarType>& a, std::pair<std::string, ScalarType>& b) {
              return std::abs(a.second) > std::abs(b.second);
            });

  maquis::cout << std::setw(log(parms_["srcas_numSamples"]) + 2) << "index" << std::setw(parms_["L"] * 2 + 2) << "ONV "
               << std::setw(20) << "Coeff." << std::setw(20) << "Completeness"
               << "\n";
  ScalarType tmpCompleteness = 0.0;
  int count = 0;
  for (const auto& i : onvValuePairs) {
    count++;
    tmpCompleteness += std::pow(i.second, 2);
    maquis::cout << std::setw(log(parms_["srcas_numSamples"]) + 2) << count << std::setw(parms_["L"] * 2 + 2) << i.first
                 << std::setw(20) << std::setprecision(14) << std::fixed << i.second << std::setw(20)
                 << std::setprecision(14) << std::fixed << tmpCompleteness << std::endl;
  }
}

template class BaseSRCAS<double, Determinant>;
template class BaseSRCAS<std::complex<double>, Determinant>;
template class BaseSRCAS<double, VibONV>;
template class BaseSRCAS<std::complex<double>, VibONV>;

} // namespace srcas
} // namespace maquis
