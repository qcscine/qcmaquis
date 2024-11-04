
#include "determinant.h"
#include <algorithm>
#include <stdexcept>

namespace maquis {
namespace srcas {

Determinant::Determinant(const std::vector<int>& onvVec) : ONV(onvVec) {
  int index = 0;
  for (const auto& orb : onvVec) {
    switch (orb) {
      case 4:
        this->alpha_occupied_.push_back(index);
        this->beta_occupied_.push_back(index);
        break;
      case 3:
        this->alpha_occupied_.push_back(index);
        this->beta_unoccupied_.push_back(index);
        break;
      case 2:
        this->alpha_unoccupied_.push_back(index);
        this->beta_occupied_.push_back(index);
        break;
      case 1:
        this->alpha_unoccupied_.push_back(index);
        this->beta_unoccupied_.push_back(index);
        break;
      default:
        throw std::runtime_error("Found determinant with unphysical occupation");
    }
    index++;
  }
}

std::vector<int> Determinant::vector() const {
  std::vector<int> detVec(nIndices(), 1);
  for (const auto& alpha : alpha_occupied_) {
    detVec[alpha] += 2;
  }
  for (const auto& beta : beta_occupied_) {
    detVec[beta] += 1;
  }
  return detVec;
}

void Determinant::excite_electron(int from, int to, Spin spin) {
  if (spin == Spin::alpha) {
    auto alpha_occ = std::find(alpha_occupied_.begin(), alpha_occupied_.end(), from);
    auto alpha_unocc = std::find(alpha_unoccupied_.begin(), alpha_unoccupied_.end(), to);
    // assert(alpha_occ);
    // assert(alpha_unocc);
    auto alpha_occ_index = alpha_occ - alpha_occupied_.begin();
    auto alpha_unocc_index = alpha_unocc - alpha_unoccupied_.begin();
    std::swap(this->alpha_occupied_[alpha_occ_index], this->alpha_unoccupied_[alpha_unocc_index]);
    return;
  }
  auto beta_occ = std::find(beta_occupied_.begin(), beta_occupied_.end(), from);
  auto beta_unocc = std::find(beta_unoccupied_.begin(), beta_unoccupied_.end(), to);
  // assert(beta_occ);
  // assert(beta_unocc);
  auto beta_occ_index = beta_occ - beta_occupied_.begin();
  auto beta_unocc_index = beta_unocc - beta_unoccupied_.begin();
  std::swap(this->beta_occupied_[beta_occ_index], this->beta_unoccupied_[beta_unocc_index]);
}

const std::vector<int>& Determinant::occupied(Spin spin) const {
  if (spin == Spin::alpha) {
    return this->alpha_occupied_;
  }
  return this->beta_occupied_;
}

const std::vector<int>& Determinant::unoccupied(Spin spin) const {
  if (spin == Spin::alpha) {
    return this->alpha_unoccupied_;
  }
  return this->beta_unoccupied_;
}

void Determinant::flipSpin() {
  // spin flips make only sense for closed shell
  if (this->alpha_occupied_.size() != this->beta_occupied_.size()) {
    return;
  }
  this->alpha_occupied_.swap(this->beta_occupied_);
  this->alpha_unoccupied_.swap(this->beta_unoccupied_);
}

bool Determinant::operator==(const Determinant& rhs) const {
  auto other = dynamic_cast<Determinant const&>(rhs);
  std::set<int> this_alpha(this->alpha_occupied_.begin(), this->alpha_occupied_.end());
  std::set<int> this_beta(this->beta_occupied_.begin(), this->beta_occupied_.end());
  std::set<int> other_alpha(other.alpha_occupied_.begin(), other.alpha_occupied_.end());
  std::set<int> other_beta(other.beta_occupied_.begin(), other.beta_occupied_.end());
  return this_alpha == other_alpha && this_beta == other_beta;
  // return this->alpha_occupied_ == other.alpha_occupied_ && this->beta_occupied_ == other.beta_occupied_;
  // It should be enough to compare both occupied vectors, as the determinant would violate particle number consistency
  // && this->alpha_unoccupied_ == rhs.alpha_unoccupied_
  // && this->beta_unoccupied_ == rhs.beta_unoccupied_;
}

} // namespace srcas
} // namespace maquis
