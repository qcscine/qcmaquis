/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

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
        maquis::cout << "Usage: mps2ci_nMode <input file> " << std::endl;
        exit(1);
    }
    maquis::cout.precision(10);
    DmrgOptions opt(argc, argv);
    if (opt.valid) {
        // Creates the simulation object
        maquis::cout.precision(10);
        if(!(opt.parms["MODEL"] == "nmode") && !(opt.parms["MODEL"] == "watson"))
            throw std::runtime_error("This app supports only vibrational Hamiltonians");
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
