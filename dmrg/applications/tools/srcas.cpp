/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2023 Reiher Group, ETH Zurich
 *               2023- by Nina Glaser <nglaser@phys.chem.ethz.ch>
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

#include "dmrg/tools/srcas_utilities.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/symmetry_factory.h"
#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"

/**
 * @brief Application that extracts the CI coefficients associated with a given MPS
 * 
 * This applications takes as input a DMRG input file, looks for the chkp file defined
 * in that input file, loads the corresponding MPS, and performs a stochastic sampling
 * of the active space to determine the CI expansion coefficients.
 */

int main(int argc, char ** argv)
{
    // Check coherence in input
    if (argc != 2) {
        maquis::cout << "Usage: srcas <input file> " << std::endl;
        exit(1);
    }
    DmrgOptions opt(argc, argv);
    if (opt.valid) {
        if(!(opt.parms["MODEL"] == "nmode") && !(opt.parms["MODEL"] == "watson") && !(opt.parms["MODEL"] == "quantum_chemistry"))
            throw std::runtime_error("The SRCAS supports only vibrational and electronic Hamiltonians so far");
        maquis::cout.precision(10);
        maquis::cout << "---------------------- SRCAS ----------------------" << std::endl << std::endl;
        // Creates the simulation object either with real or complex coefficients
        if (opt.parms["COMPLEX"]) {
            using ScalarType = std::complex<double>;
            using InterfaceType = maquis::DMRGInterface<ScalarType>;
            std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(opt.parms);
            SRCAS<ScalarType> srcas(opt.parms, interface);
            srcas.printSRCASSettings();
            srcas.run();
            srcas.printResults();
        } else {
            using ScalarType = double;
            using InterfaceType = maquis::DMRGInterface<ScalarType>;
            std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(opt.parms);
            SRCAS<ScalarType> srcas(opt.parms, interface);
            srcas.printSRCASSettings();
            srcas.run();
            srcas.printResults();
        }
    }
    else {
        throw std::runtime_error("Parameters in inputfile corrupted");
    }
    maquis::cout << std::endl;
    return 0;
}

        