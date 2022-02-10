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

#ifndef TEST_BENZENE_FIXTURE_H
#define TEST_BENZENE_FIXTURE_H

#include "maquis_dmrg.h"
#include "dmrg/block_matrix/symmetry.h"

/** @brief Fixture class constaining the input parameters for a 
 *         transcorrelated Hamiltonian */
struct TranscorrelatedFixture
{
    /** @brief Constructor for the fixture class */
    TranscorrelatedFixture() {
        // 2x2 Real-Space Fermi-Hubbard Hamiltonian.
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("L", 4);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("width_FermiHubbard", 2);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("height_FermiHubbard", 2);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("max_bond_dimension", 50);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("U_FermiHubbard", 4.);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("symmetry", "2u1");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("site_types", "0,0,0,0");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("MODEL", "fermi_hubbard_real");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("u1_total_charge1", 2);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("u1_total_charge2", 1);
        // 2x2 Fermi-Hubbard Hamiltonian.
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("L", 4);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("width_FermiHubbard", 2);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("height_FermiHubbard", 2);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("max_bond_dimension", 50);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("U_FermiHubbard", 4.);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("symmetry", "2u1");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("site_types", "0,0,0,0");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("MODEL", "fermi_hubbard_momentum");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("u1_total_charge1", 2);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("u1_total_charge2", 1);
    }

    // Class members
    DmrgParameters parameters2x2_RealSpace_U4_2Alpha1Beta, parameters2x2_MomentumSpace_U4_2Alpha1Beta;
};

#endif