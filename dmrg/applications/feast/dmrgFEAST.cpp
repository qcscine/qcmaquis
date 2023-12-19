/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "utils/io.hpp" // has to be first include because of impi
#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>
#include <sys/time.h>
#include "utils/data_collector.hpp"
#include "utils/timings.h"
#include "maquis_dmrg.h"
#include "dmrg/utils/DmrgOptions.h"


int main(int argc, char ** argv)
{
  std::cout << "  SCINE QCMaquis " << std::endl
            << "  Quantum Chemical Density Matrix Renormalization group" << std::endl
            << "  available from https://scine.ethz.ch/download/qcmaquis" << std::endl
            << "  based on the ALPS MPS codes from http://alps.comp-phys.org" << std::endl
            << "  copyright (c) 2015-2018 Department of Chemistry and Applied Biosciences, ETH Zurich" << std::endl
            << "  copyright (c) 2012-2016 by Sebastian Keller" << std::endl
            << "  copyright (c) 2016-2021 by Alberto Baiardi, Leon Freitag," << std::endl
            << "  Stefan Knecht, Yingjin Ma" << std::endl
            << "  For details on DMRG[FEAST] see the publication:" << std::endl
            << "  A. Baiardi, A.~K.~Kelemen, M.~Reiher, J. Chem. Theory Comput. 18, 415 (2022)" << std::endl;
  DmrgOptions opt(argc, argv);
  if (opt.valid) {
    DCOLLECTOR_SET_SIZE(gemm_collector, opt.parms["max_bond_dimension"]+1)
    DCOLLECTOR_SET_SIZE(svd_collector, opt.parms["max_bond_dimension"]+1)
    Timer sim("FEAST Simulation");
    sim.begin();
    maquis::DMRGInterface<std::complex<double>> interface(opt.parms);
    interface.runFEAST();
    sim.end();
    DCOLLECTOR_SAVE_TO_FILE(gemm_collector, "collectors.h5", "/results")
    DCOLLECTOR_SAVE_TO_FILE(svd_collector, "collectors.h5", "/results")
  }
}
