/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include <cmath>
#include <iterator>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <sys/stat.h>
#include <vector>

#include "maquis_dmrg.h"
#include "dmrg/sim/symmetry_factory.h"
#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"

/**
 * @brief Application that extracts the CI coefficients associated with a given MPS
 * 
 * This applications takes as input a DMRG input file, looks for the chkp file defined
 * in that input file, loads the corresponding MPS, and calculates the overlap with
 * the Slater determinants that are listed in a given input file, which is specified
 * in the DMRG input file.
 */

int main(int argc, char ** argv)
{
    // Check coherence in input
    if (argc != 2) {
        maquis::cout << "Usage: mps2ci_vib <input file> " << std::endl;
        exit(1);
    }
    DmrgOptions opt(argc, argv);
    if (opt.valid) {
        if(!(opt.parms["MODEL"] == "nmode") && !(opt.parms["MODEL"] == "watson"))
            throw std::runtime_error("This app supports only vibrational Hamiltonians");
        maquis::cout.precision(10);
        // Creates the simulation object
        maquis::DMRGInterface<double> interface(opt.parms);
        // Opens the determinant file and loops over it
        std::string nameOfDetFile = opt.parms["determinant_file"];
        double threshold = opt.parms["determinant_threshold"];
        std::ifstream is(nameOfDetFile);
        std::string str;
        while (getline(is, str)) {
            auto overlap = interface.getCICoefficient(str);
            if (std::abs(overlap) > threshold)
                std::cout << "CI coefficient of " << str << " : " << overlap << std::endl;
        }
    }
    else {
        throw std::runtime_error("Parameters file corrupted");
    }
    maquis::cout << std::endl;
    return 0;
}
