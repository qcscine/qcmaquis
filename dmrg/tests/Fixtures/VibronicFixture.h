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

#ifndef TEST_VIBRONIC_FIXTURE_H
#define TEST_VIBRONIC_FIXTURE_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

/**
 * @brief Fixture class for the test of the Vibronic DMRG
 */
struct VibronicFixture
{
    /** @brief Constructor for the fixture class */
    VibronicFixture() {
        // Excitonic models with 6 monomers, 2 electronic states, 10 modes --> L=72 
        parametersExcitonic.set("L", 72);
        parametersExcitonic.set("symmetry", "u1");
        parametersExcitonic.set("LATTICE", "vibronic lattice");
        parametersExcitonic.set("MODEL", "excitonic");
        parametersExcitonic.set("Nmax", 6);
        parametersExcitonic.set("vibronic_nstates", 2);
        parametersExcitonic.set("vibronic_nmodes", 10);
        parametersExcitonic.set("n_excitons", 6);
        // Vibronic models for S1/S2 states of pyrazine --> L=26
        parametersVibronic.set("L", 26);
        parametersVibronic.set("symmetry", "u1");
        parametersVibronic.set("LATTICE", "vibronic lattice");
        parametersVibronic.set("MODEL", "vibronic");
        parametersVibronic.set("Nmax", 6);
        parametersVibronic.set("vibronic_nstates", 2);
        parametersVibronic.set("vibronic_nmodes", 24);
        // "Fake" vibronic Hamiltonian with only a single state, with an Harmonic PES.
        parametersFakeVibronic.set("L", 4);
        parametersFakeVibronic.set("symmetry", "u1");
        parametersFakeVibronic.set("LATTICE", "vibronic lattice");
        parametersFakeVibronic.set("MODEL", "vibronic");
        parametersFakeVibronic.set("Nmax", 6);
        parametersFakeVibronic.set("vibronic_nstates", 1);
        parametersFakeVibronic.set("vibronic_nmodes", 3);
        parametersFakeVibronic.set("integral_file", "integral_file_VibronicFile");
        //
        integralFileFakeVibronic.open("integral_file_VibronicFile");
        integralFileFakeVibronic << "EL_ST 0 0 " << std::endl;
        integralFileFakeVibronic << " 1.0000   1   1  " << std::endl;
        integralFileFakeVibronic << "-1.0000  -1  -1  " << std::endl;
        integralFileFakeVibronic << " 2.0000   2   2  " << std::endl;
        integralFileFakeVibronic << "-2.0000  -2  -2  " << std::endl;
        integralFileFakeVibronic << " 3.0000   3   3  " << std::endl;
        integralFileFakeVibronic << "-3.0000  -3  -3  " << std::endl;
        integralFileFakeVibronic.close();
    }

    /** @brief Class destructor */
    ~VibronicFixture() {
        std::remove("integral_file_VibronicFile");
    }

    // Class members
    DmrgParameters parametersVibronic, parametersExcitonic, parametersFakeVibronic;
    std::ofstream integralFileFakeVibronic;
};

#endif