#ifndef MAQUIS_SRCAS_VIB_SRCAS_H
#define MAQUIS_SRCAS_VIB_SRCAS_H

#include "base_srcas.h"
#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"
#include "vib_onv.h"
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <memory>
#include <vector>

namespace maquis {
namespace srcas {

/**
 * @brief SRCAS class responsible for the sampling of an MPS
 *
 * @tparam ScalarType double or complex
 */
template<typename ScalarType>
class VibSRCAS : public virtual BaseSRCAS<ScalarType, VibONV> {
  using InterfaceType = maquis::DMRGInterface<ScalarType>;

 public:
  VibSRCAS(DmrgParameters& parameters, std::shared_ptr<InterfaceType> interface);
  /** @brief generate a new ONV form current queen */
  VibONV generateNewONV_() override;
  VibONV generateSymmetricDeterminant_(const VibONV& onv) const override;
  void setQueenFromString_(const std::string& queen) override;
  /** @brief get a random virtual orbital for det
   *
   * @param det vector representation of a Determinant
   **/
  // bool symmetriesFulfilled_(const std::vector<int>& det) const override;
  bool validONV_(const VibONV& onv) const override;

 private:
  /** @brief modes or electrons **/
  std::vector<int> onvSpace_;
};

} // namespace srcas
} // namespace maquis

#endif
