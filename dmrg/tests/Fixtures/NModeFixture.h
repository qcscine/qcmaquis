/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef TEST_NMODE_FIXTURE_H
#define TEST_NMODE_FIXTURE_H

#include "dmrg/utils/DmrgParameters.h"

/**
 * @brief Fixture class for the test of the n-mode vibrational DMRG code.
 */
struct NModeFixture
{
    /** @brief Constructor for the fixture class */
    NModeFixture() {
        // == PARAMETERS FOR DUMMY CALCULATIONS ==
        parametersTwoMode.set("L", 24);
        parametersTwoMode.set("nmode_num_modes", 2);
        parametersTwoMode.set("nmode_max_coupling", 2);
        parametersTwoMode.set("nmode_num_basis", "12,12");
        //
        parametersFourMode.set("L", 18);
        parametersFourMode.set("nmode_num_modes", 4);
        parametersFourMode.set("nmode_max_coupling", 3);
        parametersFourMode.set("nmode_num_basis", "3,4,5,6");
    }
    // Class members
    DmrgParameters parametersTwoMode, parametersFourMode;
};

#endif