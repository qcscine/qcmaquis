#ifndef MAQUIS_SRCAS_BASE_SRCAS_H
#define MAQUIS_SRCAS_BASE_SRCAS_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <memory>
#include <string>

namespace maquis {
namespace srcas {

/** @brief Virtual SRCAS class */
template<typename ScalarType, class T>
class BaseSRCAS {
  using InterfaceType = maquis::DMRGInterface<ScalarType>;

 public:
  virtual ~BaseSRCAS() = default;
  BaseSRCAS(DmrgParameters& parameters, std::shared_ptr<InterfaceType> interface);
  /** @brief Run SRCAS sampling */
  void restart();
  /** @brief Run SRCAS sampling */
  void run();
  /** @brief Run SRCAS sampling with extra determinants */
  void run(const std::vector<std::string>& extra_onvs);
  /**
   * @brief get last determinant queen
   *
   * @return vector representation of current queen
   **/
  std::vector<int> currentQueen() const;
  /**
   * @brief getter map with all sampled determinants above the threshold
   *
   * @retrun map with vector representation as key and correpsonding value
   **/
  const std::map<std::string, ScalarType>& sampledTable() const;
  /**
   * @brief get achieved completeness
   *
   * @return double the completeness
   */
  double completeness() const;
  /** @brief print general settings used in srcas and from printSettingImpl_ */
  void printSettings() const;
  /** @brief print general settings used in srcas and from printSettingImpl_ */
  void printResults() const;

 protected:
  /** @brief sample user specified determinants */
  void evaluateSpecificONVs_(const std::vector<std::string>& extra_onvs);

  virtual void setQueenFromString_(const std::string& queen) = 0;

  /** @brief evaluate the initial queen */
  void evaluateInitial_();
  /** @brief main sampling*/
  void sample_();

  /** @brief generate a new determinant form current queen */
  virtual T generateNewONV_() = 0;
  virtual T generateSymmetricDeterminant_(const T& onv) const = 0;
  /** @brief override to add settings to print in srcas settings */
  void printSettingImpl_() const;
  /** @brief get a random virtual orbital for det
   *
   * @param det vector representation of a Determinant
   **/
  virtual bool validONV_(const T& onv) const = 0;
  /** @brief boost random number generator **/
  boost::mt19937 generator_;
  /** @brief boost uniform distribution **/
  boost::uniform_real<> uniformDist_;
  /** @brief boost geometric distribution **/
  boost::geometric_distribution<double> geomDist_;
  /** @brief boost uniform distribution generator **/
  mutable boost::variate_generator<boost::mt19937&, boost::uniform_real<double>> uniformRandomNumber_;
  /** @brief boost geometric distribution generator **/
  mutable boost::variate_generator<boost::mt19937&, boost::geometric_distribution<double>> geometricRandomNumber_;
  /** @brief all DMRG parameters **/
  DmrgParameters& parms_;
  /** @brief DMRG interface **/
  std::shared_ptr<InterfaceType> interface_;
  /** @brief current completeness **/
  double completeness_;
  /** @brief all sampled coeffs above threshold **/
  std::map<std::string, ScalarType> hashTable_;
  /** @brief current queen **/
  T queen_;
  typename std::map<std::string, ScalarType>::iterator iter_;

  // useful for restart
  /** @brief indicate the start of the sampling */
  std::string startString_ = "----- Starting SRCAS -----";
  /** @brief indicate the end of the sampling */
  std::string endString_ = "----- End SRCAS -----";
};

} // namespace srcas
} // namespace maquis

#endif
