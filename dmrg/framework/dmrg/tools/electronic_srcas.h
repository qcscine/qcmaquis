#ifndef MAQUIS_SRCAS_ELECTRONIC_SRCAS_H
#define MAQUIS_SRCAS_ELECTRONIC_SRCAS_H

#include "base_srcas.h"
#include "determinant.h"
#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"
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
class ElectronicSRCAS : public virtual BaseSRCAS<ScalarType, Determinant> {
  using InterfaceType = maquis::DMRGInterface<ScalarType>;

 public:
  ElectronicSRCAS(DmrgParameters& parameters, std::shared_ptr<InterfaceType> interface);
  /** @brief generate a new determinant form current queen */
  Determinant generateNewONV_() override;
  /** @brief get a random occupied orbital for det
   *
   * @param det vector representation of a Determinant
   **/
  int getRandomOccSpinOrb_(const Determinant& det, Spin spin) const;
  /** @brief get a random virtual orbital for det
   *
   * @param det vector representation of a Determinant
   **/
  int getRandomUnoccSpinOrb_(const Determinant& det, Spin spin) const;
  /** @brief get a random virtual orbital for det
   *
   * @param det vector representation of a Determinant
   **/
  // bool symmetriesFulfilled_(const std::vector<int>& det) const override;
  bool validONV_(const Determinant& det) const override;

 private:
  /** @brief modes or electrons **/
  double fractional_beta_orbs_;
};

} // namespace srcas
} // namespace maquis

#endif
