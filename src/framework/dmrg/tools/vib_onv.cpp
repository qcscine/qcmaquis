#include "vib_onv.h"

namespace maquis {
namespace srcas {

VibONV::VibONV(const std::vector<int>& onvVec
)  ///{ // , const std::vector<int>& onvSpace)
    : ONV(onvVec), onv_(onvVec) {}

std::vector<int> VibONV::vector() const { return onv_; }

bool VibONV::operator==(const VibONV& rhs) const {
  return this->onv_ == rhs.onv_;
  // It should be enough to compare both occupied vectors, as the determinant
  // would violate particle number consistency
  // && this->alpha_unoccupied_ == rhs.alpha_unoccupied_
  // && this->beta_unoccupied_ == rhs.beta_unoccupied_;
}

}  // namespace srcas
}  // namespace maquis
