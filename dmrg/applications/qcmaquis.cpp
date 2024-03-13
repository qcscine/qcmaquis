/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#include <sys/stat.h>
#include <sys/time.h>

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "dmrg/utils/DmrgOptions.h"
#include "maquis_dmrg.h"
#include "utils/data_collector.hpp"
#include "utils/timings.h"

// @brief Checks that request simulation has been enabled at compile time
void checkEnabledSimulationType(const std::string& sim_type) {
#ifndef DMRG_TD
  if (sim_type == "evolve") {
    std::cerr << "Time-dependent DMRG not available. Please recompile "
                 "with -DBUILD_DMRG_EVOLVE=ON\n";
    exit(1);
  }
#endif
#ifndef DMRG_FEAST
  if (sim_type == "feast") {
    std::cerr << "FEAST not available. Please recompile with "
                 "-DBUILD_DMRG_FEAST=ON\n";
    exit(1);
  }
#endif
#ifndef DMRG_TRANSCORRELATED
  if (sim_type == "transcorrelated") {
    std::cerr << "Transcorrelated DMRG not available. Please recompile "
                 "with -DBUILD_TRANSCORRELATED_DMRG=ON\n";
    exit(1);
  }
#endif
}

// @brief Checks that the requested simulation type is valid
void checkSimulationType(const std::string& sim_type) {
  std::vector<std::string> valid_types = {
      "optimize", "evolve", "ipi", "feast", "transcorrelated"
  };
  if (std::find(valid_types.begin(), valid_types.end(), sim_type) ==
      valid_types.end()) {
    std::cerr << "Unknown simulation type: \"" << sim_type
              << "\". Valid options are: optimize, evolve, ipi, feast, "
                 "transcorrelated\n";
    exit(1);
  }
  checkEnabledSimulationType(sim_type);
}

int main(int argc, char** argv) {
  std::cout << "  SCINE QCMaquis \n"
            << "  Quantum Chemical Density Matrix Renormalization group\n"
            << "  available from https://scine.ethz.ch/download/qcmaquis\n"
            << "  based on the ALPS MPS codes from http://alps.comp-phys.org/\n"
            << "  copyright (c) 2015-2018 Department of Chemistry and Applied "
               "Biosciences, ETH Zurich\n"
            << "  copyright (c) 2012-2016 by Sebastian Keller\n"
            << "  copyright (c) 2016-2018 by Alberto Baiardi, Leon Freitag, \n"
            << "  Stefan Knecht, Yingjin Ma \n"
            << "  for details see the publication: \n"
            << "  S. Keller et al., J. Chem. Phys. 143, 244118 (2015)\n"
            << std::endl;

  DmrgOptions opt(argc, argv);

  if (opt.valid) {
    if (opt.parms.is_set("simulation_type")) {
      std::string sim_type = opt.parms["simulation_type"];
      checkSimulationType(sim_type);

      DCOLLECTOR_SET_SIZE(gemm_collector, opt.parms["max_bond_dimension"] + 1)
      DCOLLECTOR_SET_SIZE(svd_collector, opt.parms["max_bond_dimension"] + 1)

      Timer sim_timer(
          "\n*********************************************************"
          "\nQCMAQUIS " +
          sim_type + " simulation"
      );
      sim_timer.begin();

      if (!opt.parms["COMPLEX"]) {
        maquis::DMRGInterface<double> interface(opt.parms);
        if (opt.parms["fiedler"]) {
          opt.parms["orbital_order"] =
              interface.fiedler_order(1, std::vector<std::vector<int>>{}, "");
        }
        interface.run(sim_type);
      } else {
        maquis::DMRGInterface<std::complex<double>> interface(opt.parms);
        if (opt.parms["fiedler"]) {
          opt.parms["orbital_order"] =
              interface.fiedler_order(1, std::vector<std::vector<int>>{}, "");
        }
        interface.run(sim_type);
      }

      sim_timer.end();
    }
  }
}
