/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#include <sys/stat.h>
#include <sys/time.h>

#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>

#include "dmrg/utils/DmrgOptions.h"
#include "maquis_dmrg.h"
#include "utils/data_collector.hpp"
#include "utils/timings.h"

int main(int argc, char **argv) {
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
      DCOLLECTOR_SET_SIZE(gemm_collector, opt.parms["max_bond_dimension"] + 1)
      DCOLLECTOR_SET_SIZE(svd_collector, opt.parms["max_bond_dimension"] + 1)
      Timer sim(sim_type);
      sim.begin();
      if (sim_type == "optimization") {
        // Here we must explicitly distinguish all cases.
        if (!opt.parms["COMPLEX"]) {
          maquis::DMRGInterface<double> interface(opt.parms);
          interface.optimize();
        } else {
          maquis::DMRGInterface<std::complex<double>> interface(opt.parms);
          interface.optimize();
        }
      } else if (sim_type == "measure") {
        if (!opt.parms["COMPLEX"]) {
          maquis::DMRGInterface<double> interface(opt.parms);
          interface.run_measure();
        } else {
          maquis::DMRGInterface<std::complex<double>> interface(opt.parms);
          interface.run_measure();
        }
      } else if (sim_type == "inverse_power_iteration") {
        if (!opt.parms["COMPLEX"]) {
          maquis::DMRGInterface<double> interface(opt.parms);
          interface.runInversePowerIteration();
        } else {
          maquis::DMRGInterface<std::complex<double>> interface(opt.parms);
          interface.runInversePowerIteration();
        }
      } else if (sim_type == "time_dep") {
#ifdef DMRG_TD
        maquis::DMRGInterface<std::complex<double>> interface(opt.parms);
        interface.evolve();
#else
        std::cerr << "Time-dependent DMRG not available. Please recompile "
                     "with -DBUILD_DMRG_EVOLVE=ON\n";
        exit(1);
#endif
      } else if (sim_type == "feast") {
#ifdef DMRG_FEAST
        maquis::DMRGInterface<std::complex<double>> interface(opt.parms);
        interface.runFEAST();
#else
        std::cerr << "FEAST not available. Please recompile with "
                     "-DBUILD_DMRG_FEAST=ON\n";
        exit(1);
#endif
      } else if (sim_type == "transcorrelation") {
#ifdef DMRG_TRANSCORRELATED
        maquis::DMRGInterface<double> interface(opt.parms);
        interface.runTranscorrelated();
#else
        std::cerr << "Transcorrelated DMRG not available. Please recompile "
                     "with -DBUILD_TRANSCORRELATED_DMRG=ON\n";
        exit(1);
#endif
      } else {
        std::cerr << "Unknown simulation type: " << sim_type
                  << ". Valid options are: optimization, measure, "
                     "inverse-power-iteration, evolve, feast\n";
        exit(1);
      }
      sim.end();
      DCOLLECTOR_SAVE_TO_FILE(gemm_collector, "collectors.h5", "/results")
      DCOLLECTOR_SAVE_TO_FILE(svd_collector, "collectors.h5", "/results")
    } else {
      std::cerr << "Simulation type not set. Valid options are: optimization, "
                   "measure, "
                   "inverse-power-iteration, evolve, feast\n";
      exit(1);
    }
  }
}
