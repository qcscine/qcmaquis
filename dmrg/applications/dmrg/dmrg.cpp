/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
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
    std::cout << "  SCINE QCMaquis \n"
              << "  Quantum Chemical Density Matrix Renormalization group\n"
              << "  available from https://scine.ethz.ch/download/qcmaquis\n"
              << "  based on the ALPS MPS codes from http://alps.comp-phys.org/\n"
              << "  copyright (c) 2015-2018 Laboratory of Physical Chemistry, ETH Zurich\n"
              << "  copyright (c) 2012-2016 by Sebastian Keller\n"
              << "  copyright (c) 2016-2018 by Alberto Baiardi, Leon Freitag, \n"
              << "  Stefan Knecht, Yingjin Ma \n"
              << "  for details see the publication: \n"
              << "  S. Keller et al., J. Chem. Phys. 143, 244118 (2015)\n"
              << std::endl;

    DmrgOptions opt(argc, argv);
    if (opt.valid) {
        maquis::cout.precision(10);
        DCOLLECTOR_SET_SIZE(gemm_collector, opt.parms["max_bond_dimension"]+1)
        DCOLLECTOR_SET_SIZE(svd_collector, opt.parms["max_bond_dimension"]+1)
        timeval now, then, snow, sthen;
        gettimeofday(&now, NULL);
        // Here we must explicitly distinguish all cases.
        if (!opt.parms["COMPLEX"]) {
            maquis::DMRGInterface<double> interface(opt.parms);
            interface.optimize();
        }
        else {
            maquis::DMRGInterface<std::complex<double> > interface(opt.parms);
            interface.optimize();
        }
        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
        DCOLLECTOR_SAVE_TO_FILE(gemm_collector, "collectors.h5", "/results")
        DCOLLECTOR_SAVE_TO_FILE(svd_collector, "collectors.h5", "/results")
        maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
    }
}

