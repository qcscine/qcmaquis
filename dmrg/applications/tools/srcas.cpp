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
#include <regex>
#include <string>

namespace detail {
template<class T>
void runSRCAS(T& srcas, const std::vector<std::string>& onvs, bool restart = false) {
  srcas.printSettings();
  if (restart) {
    srcas.restart();
  }
  srcas.run(onvs);
  srcas.printResults();
}

template<class T>
void runSRCAS(T& srcas, bool restart = false) {
  srcas.printSettings();
  if (restart) {
    srcas.restart();
  }
  srcas.run();
  srcas.printResults();
}

std::vector<std::string> readONVsFromFile(const std::string& filename) {
  std::ifstream config_file;
  config_file.open(filename.c_str());
  if (!config_file.is_open()) {
    maquis::cout << "Can not open file: " << filename << std::endl;
    exit(1);
  }

  // Match only valid determinants like:
  // 4,3,2,1
  // 4,3,2,1
  // 4,3,2,1
  // But does not match
  // 4,3,2,5
  // 4,3,2,1,
  // Does also not check number of orbitals
  std::regex det_reg("^(([1234],)+[1234]{1})$");
  std::string line;
  std::smatch det_match;
  std::vector<std::string> onvs;
  int line_num = 0;

  while (std::getline(config_file, line)) {
    line_num++;
    // remove leading and trailing whitespaces
    line = std::regex_replace(line, std::regex("^ +| +$"), "$1");

    if (std::regex_match(line, det_match, det_reg)) {
      onvs.push_back(line);
    }
    else {
      maquis::cout << "Found invalid ONV in line " << line_num << ": " << line << std::endl;
    }
  }

  // remove duplicate onvs
  onvs.erase(std::unique(onvs.begin(), onvs.end()), onvs.end());

  for (auto& det : onvs) {
    maquis::cout << "read onv: <" << det << ">" << std::endl;
  }

  return onvs;
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
  maquis::cout << "\n---------------------------------------------------\n";
  maquis::cout << "-                                                 -\n";
  maquis::cout << "-                      SRCAS                      -\n";
  maquis::cout << "-                                                 -\n";
  maquis::cout << "---------------------------------------------------\n\n";
  bool restart = false;
  if (opt.parms.is_set("srcas_restart")) {
    restart = true;
    // maquis::cout << "NOT IMPLEMENTED YET!!!!" << std::endl;
    maquis::cout << "Restart from " << opt.parms["srcas_restart"] << std::endl;
  }

  if (opt.parms["COMPLEX"]) {
    using ScalarType = std::complex<double>;
    using InterfaceType = maquis::DMRGInterface<ScalarType>;
    std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(opt.parms);

    if ((opt.parms["MODEL"] == "quantum_chemistry")) {
      maquis::srcas::ElectronicSRCAS<ScalarType> srcas(opt.parms, interface);
      if (opt.parms.is_set("srcas_detfile")) {
        maquis::cout << "Reading determinants from file " << opt.parms["srcas_detfile"] << std::endl;
        auto onvs = detail::readONVsFromFile(opt.parms["srcas_detfile"]);
        detail::runSRCAS(srcas, onvs, restart);
      }
      else {
        detail::runSRCAS(srcas, restart);
      }
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
    if (opt.parms.is_set("srcas_detfile")) {
      auto onvs = detail::readONVsFromFile(opt.parms["srcas_detfile"]);
      detail::runSRCAS(srcas, onvs, restart);
    }
    else {
      detail::runSRCAS(srcas, restart);
    }
  }

  else {
    maquis::srcas::VibSRCAS<ScalarType> srcas(opt.parms, interface);
    detail::runSRCAS(srcas);
  }

  return 0;
}
