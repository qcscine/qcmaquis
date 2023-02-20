/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "utils/io.hpp" // has to be first include because of impi

#include "maquis_dmrg.h"
#include "dmrg/utils/DmrgOptions.h"

#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

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

        timeval now, then, snow, sthen;
        gettimeofday(&now, NULL);

        if (!opt.parms["COMPLEX"])
        {
            maquis::DMRGInterface<double> interface(opt.parms);
            interface.run_measure();
        }
        else
        {
            maquis::DMRGInterface<std::complex<double> > interface(opt.parms);
            interface.run_measure();
        }

        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);

        maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
    }
}
