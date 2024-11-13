#ifndef MAQUIS_SRCAS_VIB_ONV_H
#define MAQUIS_SRCAS_VIB_ONV_H

#include "onv.h"
#include <vector>

namespace maquis {
namespace srcas {

/** @brief Base Occupation Number Vector class. */
class VibONV : public virtual ONV {
 public:
  /** @brief do not use this class */
  VibONV() = default;
  /** @brief do not use this class */
  VibONV(const std::vector<int>& onvVec
  );  // , const std::vector<int>& onvSpace); //  : nIndices_(onvVec.size()) {
  /** @brief return the qcmaquis compatible vector repr. of the ONV. */
  std::vector<int> vector() const override;
  /** @brief compare two onvs */
  bool operator==(const VibONV& rhs) const;

 private:
  std::vector<int> onv_;
  // std::vector<int> onvSpace_;
};

}  // namespace srcas
}  // namespace maquis

#endif
