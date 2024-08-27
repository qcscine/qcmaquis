#include "electronic_srcas.h"

namespace maquis {
namespace srcas {

template<typename ScalarType>
int ElectronicSRCAS<ScalarType>::getRandomOccSpinOrb_(const Determinant& det, Spin spin) const {
  // TODO: make index choice dependent on s1 entropies if available
  // Weight entropies based on max entropy, with min entropy
  auto occupied_orbs = det.occupied(spin);
  int whichIndex = int(floor(this->uniformRandomNumber_() * occupied_orbs.size()));
  return occupied_orbs[whichIndex];
}

template<typename ScalarType>
int ElectronicSRCAS<ScalarType>::getRandomUnoccSpinOrb_(const Determinant& det, Spin spin) const {
  // TODO: make index choice dependent on s1 entropies if available
  auto unoccupied_orbs = det.unoccupied(spin);
  int whichIndex = int(floor(this->uniformRandomNumber_() * unoccupied_orbs.size()));
  return unoccupied_orbs[whichIndex];
}

template<typename ScalarType>
Determinant ElectronicSRCAS<ScalarType>::generateSymmetricDeterminant_(const Determinant& onv) const {
  auto tmpDet = onv;
  tmpDet.flipSpin();
  return tmpDet;
}

template<typename ScalarType>
Determinant ElectronicSRCAS<ScalarType>::generateNewONV_() {
  // Start from queen
  auto tmpDet = this->queen_;
  // while (tmpDet == this->queen_) {
  int nele_excited = this->geometricRandomNumber_();
  for (int i = 0; i < nele_excited; i++) {
    Spin spin = Spin::alpha;
    if (this->uniformRandomNumber_() < fractional_beta_orbs_) {
      spin = Spin::beta;
    }
    int annihilate = this->getRandomOccSpinOrb_(tmpDet, spin);
    int create = this->getRandomUnoccSpinOrb_(tmpDet, spin);
    tmpDet.excite_electron(annihilate, create, spin);
  }
  return tmpDet;
}

template<typename ScalarType>
void ElectronicSRCAS<ScalarType>::setQueenFromString_(const std::string& queen) {
  std::vector<int> queenVec;
  for (const char& i : queen) {
    if (i != ',') {
      // convert i to int
      queenVec.push_back(i - '0');
    }
  }
  this->queen_ = Determinant(queenVec);
}

template<typename ScalarType>
ElectronicSRCAS<ScalarType>::ElectronicSRCAS(DmrgParameters& parameters, std::shared_ptr<InterfaceType> interface)
  : BaseSRCAS<ScalarType, Determinant>(parameters, interface) {
  int alpha = int(this->parms_["u1_total_charge1"]);
  int beta = int(this->parms_["u1_total_charge2"]);
  fractional_beta_orbs_ = static_cast<double>(beta) / static_cast<double>(alpha + beta);

  if (this->parms_.is_set("init_basis_state")) {
    this->queen_ = Determinant(this->parms_["init_basis_state"].template as<std::vector<int>>());
    // TODO: check if det is valid
  }
  else {
    std::vector<int> tmpVec(this->parms_["L"], 1);
    for (int i = 0; i < alpha; i++) {
      tmpVec[i] += 2;
    }
    for (int i = 0; i < beta; i++) {
      tmpVec[i] += 1;
    }
    this->queen_ = Determinant(tmpVec);
  }
}

// TODO: Write this function
template<typename ScalarType>
bool ElectronicSRCAS<ScalarType>::validONV_(const Determinant& det) const {
  return true;
}

// TODO:
/* think about restaring:
 * What to choose as new queen? -> Last written determinant
 * start from last n iteration printed
 * Choose same seed, skip n iterations (shold same as seed is same)
 */

// Explicit template instantiation
template class ElectronicSRCAS<double>;
template class ElectronicSRCAS<std::complex<double>>;

} // namespace srcas
} // namespace maquis
