
#include "onv.h"

namespace maquis {
namespace srcas {

std::string ONV::string() const {
  std::vector<int> detVec = this->vector();
  std::string detString;
  for (const auto& orb : detVec) {
    detString += std::to_string(orb);
    detString += ",";
  }
  detString.pop_back();
  return detString;
}

size_t ONV::nIndices() const {
  return nIndices_;
}

bool ONV::operator==(const ONV& rhs) const {
  return this->nIndices_ == rhs.nIndices_;
}

} // namespace srcas
} // namespace maquis
