/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#include "dmrg/tools/base_srcas.h"
#include "dmrg/tools/electronic_srcas.h"
#include "dmrg/tools/vib_srcas.h"
#include "dmrg/utils/DmrgOptions.h"
#include "maquis_dmrg.h"
#include <sys/stat.h>
#include <sys/time.h>
#include <iostream>
#include <string>

namespace detail {
template<class T>
void runSRCAS(T& srcas) {
  srcas.printSettings();
  srcas.run();
  srcas.printResults();
}
} // namespace detail

/**
 * @brief Application that extracts the CI coefficients associated with a given
 * MPS
 *
 * This applications takes as input a DMRG input file, looks for the chkp file
 * defined in that input file, loads the corresponding MPS, and performs a
 * stochastic sampling of the active space to determine the CI expansion
 * coefficients.
 */
int main(int argc, char** argv) {
  // Check coherence in input
  if (argc != 2) {
    maquis::cout << "Usage: srcas <input file>\n";
    exit(1);
  }
  DmrgOptions opt(argc, argv);
  if (!opt.valid) {
    maquis::cout << "DMRG options are not valid";
    exit(1);
  }
  if (!(opt.parms["MODEL"] == "nmode") && !(opt.parms["MODEL"] == "watson") && !(opt.parms["MODEL"] == "quantum_chemistry")) {
    maquis::cout << "SRCAS is not implemented for model: " << opt.parms["MODEL"];
    exit(1);
  }

  maquis::cout.precision(10);
  maquis::cout << "---------------------- SRCAS ----------------------\n\n";

  if (opt.parms["COMPLEX"]) {
    using ScalarType = std::complex<double>;
    using InterfaceType = maquis::DMRGInterface<ScalarType>;
    std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(opt.parms);
    if ((opt.parms["MODEL"] == "quantum_chemistry")) {
      maquis::srcas::ElectronicSRCAS<ScalarType> srcas(opt.parms, interface);
      detail::runSRCAS(srcas);
    }
    else {
      maquis::srcas::VibSRCAS<ScalarType> srcas(opt.parms, interface);
      detail::runSRCAS(srcas);
    }
    return 0;
  }

  using ScalarType = double;
  using InterfaceType = maquis::DMRGInterface<ScalarType>;
  std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(opt.parms);
  if ((opt.parms["MODEL"] == "quantum_chemistry")) {
    maquis::srcas::ElectronicSRCAS<ScalarType> srcas(opt.parms, interface);
    detail::runSRCAS(srcas);
  }
  else {
    maquis::srcas::VibSRCAS<ScalarType> srcas(opt.parms, interface);
    detail::runSRCAS(srcas);
  }

  return 0;
}
