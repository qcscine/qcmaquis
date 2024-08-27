#include "vib_srcas.h"
#include <stdexcept>

namespace maquis {
namespace srcas {

template<typename ScalarType>
VibSRCAS<ScalarType>::VibSRCAS(DmrgParameters& parameters, std::shared_ptr<InterfaceType> interface)
  : BaseSRCAS<ScalarType, VibONV>(parameters, interface) {
  // std::vector<int> onvSpace;
  int nParticles = 0;
  if (this->parms_["MODEL"] == "nmode") {
    nParticles = this->parms_["nmode_num_modes"];
    onvSpace_ = this->parms_["nmode_num_basis"].template as<std::vector<int>>();
  }
  else {
    nParticles = this->parms_["L"];
    onvSpace_ = this->parms_["Nmax"].template as<std::vector<int>>();
    if (onvSpace_.size() != nParticles && onvSpace_.size() != 1) {
      throw std::runtime_error("The Nmax parameter must be either a single integer, or a vector of lenght L");
    }

    if (onvSpace_.size() != nParticles) {
      // TODO: Check if this is correct
      onvSpace_ = std::vector<int>(nParticles, std::stoi(this->parms_["Nmax"].str()));
    }
  }

  if (this->parms_.is_set("init_basis_state")) {
    this->queen_ = VibONV(this->parms_["init_basis_state"].template as<std::vector<int>>());
  }
  else {
    std::vector<int> tmpVec(nParticles, 0);
    this->queen_ = VibONV(tmpVec);
  }
  if (!this->validONV_(this->queen_.vector())) {
    std::string tmp = "Initial basis state <" + this->queen_.string() + "> is not valid";
    throw std::runtime_error(tmp);
  }
}

template<typename ScalarType>
VibONV VibSRCAS<ScalarType>::generateSymmetricDeterminant_(const VibONV& onv) const {
  throw std::runtime_error("Should always be catched");
}

template<typename ScalarType>
void VibSRCAS<ScalarType>::setQueenFromString_(const std::string& queen) {
  throw std::runtime_error("Not implemented for VibSRCAS");
}

template<typename ScalarType>
bool VibSRCAS<ScalarType>::validONV_(const VibONV& onv) const {
  auto onvVec = onv.vector();
  for (int i = 0; i < onvVec.size(); i++) {
    if (onvVec[i] >= onvSpace_[i]) {
      return false;
    }
  }
  return true;
}

template<typename ScalarType>
VibONV VibSRCAS<ScalarType>::generateNewONV_() {
  auto tmpVec = this->queen_.vector();
  while (tmpVec == this->queen_.vector()) {
    for (int i = 0; i < tmpVec.size(); i++) {
      // poisson distribution centered on the current modal
      boost::poisson_distribution<> poissonDist(tmpVec[i] + 0.5);
      boost::variate_generator<boost::mt19937&, boost::poisson_distribution<>> poissonRandomNumber(this->generator_, poissonDist);

      while (true) {
        if (this->uniformRandomNumber_() < this->parms_["srcas_samplingSpeed"]) {
          tmpVec[i] = poissonRandomNumber();
        }
        if (tmpVec[i] < onvSpace_[i]) {
          break;
        }
      }
    }
  }
  return VibONV(tmpVec);
}

template class VibSRCAS<double>;
template class VibSRCAS<std::complex<double>>;

} // namespace srcas
} // namespace maquis
