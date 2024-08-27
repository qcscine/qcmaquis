#ifndef MAQUIS_SRCAS_ONV_H
#define MAQUIS_SRCAS_ONV_H

#include <string>
#include <vector>

namespace maquis {
namespace srcas {

/** @brief Base Occupation Number Vector class. */
class ONV {
 public:
  /** @brief do not use this class */
  ONV() = default;
  /** @brief do not use this class */
  ONV(const std::vector<int>& onvVec) : nIndices_(onvVec.size()) {
  }
  /** @brief do not use this class*/
  virtual ~ONV() = default;
  /** @brief return the qcmaquis compatible vector repr. of the ONV. */
  virtual std::vector<int> vector() const = 0;
  /** @brief return the qcmaquis compatible string repr. of the ONV. */
  std::string string() const;
  /** @brief compare two onvs */
  bool operator==(const ONV& rhs) const;
  /** @brief return number of indices in the ONV*/
  size_t nIndices() const;

 private:
  /** @brief number of indices in ONV*/
  size_t nIndices_;
};

} // namespace srcas
} // namespace maquis

#endif
