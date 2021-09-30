/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef TEST_WATSON_FIXTURE_H
#define TEST_WATSON_FIXTURE_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

/**
 * @brief Fixture class for the test of the Watson Hamiltonian-based
 * vibrational DMRG code.
 */
struct WatsonFixture
{
    // Types definition
    using MaquisIntegralType = maquis::integral_map<double, chem::Hamiltonian::VibrationalCanonical>;

    /** @brief Constructor for the fixture class */
    WatsonFixture() {
        // == INPUT FILE CREATIONS ==
        integralFileEthyleneHarmonic.open("integral_file_test_Watson_Ethylene_Harmonic");
        integralFileEthyleneHarmonic << " 2.06242167E+02   1   1  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-2.06242167E+02  -1  -1  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 2.37547247E+02   2   2  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-2.37547247E+02  -2  -2  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 2.41596974E+02   3   3  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-2.41596974E+02  -3  -3  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 2.62702319E+02   4   4  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-2.62702319E+02  -4  -4  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 3.11690974E+02   5   5  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-3.11690974E+02  -5  -5  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 3.42345676E+02   6   6  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-3.42345676E+02  -6  -6  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 3.69620153E+02   7   7  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-3.69620153E+02  -7  -7  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 4.18142252E+02   8   8  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-4.18142252E+02  -8  -8  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 7.85228686E+02   9   9  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-7.85228686E+02  -9  -9  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 7.89209319E+02  10  10  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-7.89209319E+02 -10 -10  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 8.05722986E+02  11  11  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-8.05722986E+02 -11 -11  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << " 8.12177325E+02  12  12  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic << "-8.12177325E+02 -12 -12  0  0  0  0  " << std::endl;
        integralFileEthyleneHarmonic.close();
        //
        parametersEthyleneWatson.set("L", 12);
        parametersEthyleneWatson.set("symmetry", "none");
        parametersEthyleneWatson.set("LATTICE", "watson");
        parametersEthyleneWatson.set("MODEL", "watson");
        parametersEthyleneWatson.set("Nmax", 6);
        parametersEthyleneWatson.set("integral_file", "integral_file_test_Watson_Ethylene_Harmonic");
    }

    /** @brief Class destructor (removes tmp files) */
    ~WatsonFixture() {
        std::remove("integral_file_test_Watson_Ethylene_Harmonic");
    }

    // Class members
    DmrgParameters parametersEthyleneWatson;
    std::ofstream integralFileEthyleneHarmonic;
};

#endif