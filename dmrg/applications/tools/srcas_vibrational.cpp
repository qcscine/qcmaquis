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

#include "srcas_utilities.h"
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
        maquis::cout << "Usage: srcas_vib <input file> " << std::endl;
        exit(1);
    }
    DmrgOptions opt(argc, argv);
    if (opt.valid) {
        if(!(opt.parms["MODEL"] == "nmode") && !(opt.parms["MODEL"] == "watson"))
            throw std::runtime_error("This app supports only vibrational Hamiltonians");
        maquis::cout.precision(10);
        maquis::cout << "---------------------- VIBRATIONAL SRCAS ----------------------" << std::endl << std::endl;
        // Creates the simulation object
        using InterfaceType = maquis::DMRGInterface<double>;
        // EDIT NINA this should be changed at some point!
        /*
        if (opt.parms["MODEL"] == "watson")
            auto interface = std::make_shared<maquis::DMRGInterface<double, Hamiltonian::VibrationalCanonical>>(opt.parms);
        else
            auto interface = std::make_shared<maquis::DMRGInterface<double, Hamiltonian::VibrationalNMode>>(opt.parms);
        */
        std::shared_ptr<InterfaceType> interface = std::make_shared<InterfaceType>(opt.parms);
        SRCAS<InterfaceType> srcas(opt.parms, interface);
        srcas.printSRCASSettings();
        srcas.run();
        srcas.printResults();
    }
    else {
        throw std::runtime_error("Parameters in inputfile corrupted");
    }
    maquis::cout << std::endl;
    return 0;
}

        