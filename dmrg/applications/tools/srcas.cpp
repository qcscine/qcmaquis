/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#include <vector>

#include "dmrg/sim/symmetry_factory.h"
#include "dmrg/tools/srcas_utilities.h"
#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

/**
 * @brief Application that extracts the CI coefficients associated with a given
 * MPS
 *
 * This applications takes as input a DMRG input file, looks for the chkp file
 * defined in that input file, loads the corresponding MPS, and performs a
 * stochastic sampling of the active space to determine the CI expansion
 * coefficients.
 */

int main(int argc, char **argv) {
  // Check coherence in input
  if (argc != 2) {
    maquis::cout << "Usage: srcas <input file> " << std::endl;
    exit(1);
  }
  DmrgOptions opt(argc, argv);
  if (opt.valid) {
    if (!(opt.parms["MODEL"] == "nmode") && !(opt.parms["MODEL"] == "watson") &&
        !(opt.parms["MODEL"] == "quantum_chemistry"))
      throw std::runtime_error("The SRCAS supports only vibrational and "
                               "electronic Hamiltonians so far");
    maquis::cout.precision(10);
    maquis::cout << "---------------------- SRCAS ----------------------"
                 << std::endl
                 << std::endl;
    // Creates the simulation object either with real or complex coefficients
    if (opt.parms["COMPLEX"]) {
      using ScalarType = std::complex<double>;
      using InterfaceType = maquis::DMRGInterface<ScalarType>;
      std::shared_ptr<InterfaceType> interface =
          std::make_shared<InterfaceType>(opt.parms);
      SRCAS<ScalarType> srcas(opt.parms, interface);
      srcas.printSRCASSettings();
      srcas.run();
      srcas.printResults();
    } else {
      using ScalarType = double;
      using InterfaceType = maquis::DMRGInterface<ScalarType>;
      std::shared_ptr<InterfaceType> interface =
          std::make_shared<InterfaceType>(opt.parms);
      SRCAS<ScalarType> srcas(opt.parms, interface);
      srcas.printSRCASSettings();
      srcas.run();
      srcas.printResults();
    }
  } else {
    throw std::runtime_error("Parameters in inputfile corrupted");
  }
  maquis::cout << std::endl;
  return 0;
}
