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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include "utils/io.hpp"
#include <iostream>
#include "maquis_dmrg.h"

///** @brief iTD-DMRG calculation on the Fermi-Hubbard model */
//BOOST_AUTO_TEST_CASE( TestImaginaryTimeFermiHubbard )
//{
//#ifdef HAVE_TwoU1
//    // Two-site evolutions
//    DmrgParameters p;
//    p.set("time_step", 10.);
//    p.set("L", 9);
//    p.set("width_FermiHubbard", 3);
//    p.set("height_FermiHubbard", 3);
//    p.set("nsweeps", 10);
//    p.set("max_bond_dimension", 300);
//    p.set("u1_total_charge1", 4);
//    p.set("u1_total_charge2", 4);
//    p.set("symmetry", "2u1");
//    p.set("MODEL", "fermi_hubbard_real");
//    p.set("site_types", "0,0,0,0,0,0,0,0,0");
//    p.set("U_FermiHubbard", 8.);
//    p.set("propagator_maxiter", 10);
//    p.set("imaginary_time", "yes");
//    p.set("TD_backpropagation", "no");
//    p.set("simulation_type", "TD");
//    p.set("COMPLEX", 1);
//    p.set("time_units", "fs");
//    // TD
//    maquis::DMRGInterface<std::complex<double>> interfaceTD(p);
//    interfaceTD.evolve();
//    auto energyTD = std::real(interfaceTD.energy());
//    // TI
//    maquis::DMRGInterface<std::complex<double>> interfaceTI(p);
//    interfaceTI.optimize();
//    auto energyTI = std::real(interfaceTI.energy());
//    BOOST_CHECK_CLOSE(energyTD, energyTI, 1.0E-8);
//#endif // HAVE_TwoU1
//}
//
#ifdef HAVE_TwoU1
//
///** @brief iTD-DMRG calculation on the 2x2 Fermi-Hubbard model */
//BOOST_AUTO_TEST_CASE( TestImaginaryTimeFermiHubbard_2x2_Conventional )
//{
//    // Two-site evolutions
//    DmrgParameters p;
//    p.set("L", 4);
//    p.set("width_FermiHubbard", 2);
//    p.set("height_FermiHubbard", 2);
//    p.set("nsweeps", 10);
//    p.set("max_bond_dimension", 10);
//    p.set("u1_total_charge1", 2);
//    p.set("u1_total_charge2", 1);
//    p.set("symmetry", "2u1");
//    p.set("MODEL", "fermi_hubbard_real");
//    p.set("site_types", "0,0,0,0");
//    p.set("U_FermiHubbard", 4.);
//    p.set("propagator_maxiter", 10);
//    p.set("imaginary_time", "yes");
//    p.set("TD_backpropagation", "no");
//    p.set("simulation_type", "TD");
//    p.set("COMPLEX", 1);
//    p.set("time_units", "fs");
//    p.set("time_step", 10.);
//    // TD
//    maquis::DMRGInterface<std::complex<double>> interfaceTD(p);
//    interfaceTD.evolve();
//    auto energyTD = std::real(interfaceTD.energy());
//    // TI
//    maquis::DMRGInterface<std::complex<double>> interfaceTI(p);
//    interfaceTI.optimize();
//    auto energyTI = std::real(interfaceTI.energy());
//    BOOST_CHECK_CLOSE(energyTD, energyTI, 1.0E-8);
//}

/** @brief iTD-DMRG calculation on the 2x2 Fermi-Hubbard model */
BOOST_AUTO_TEST_CASE( TestImaginaryTimeFermiHubbard_2x2_Transcorrelated_1 )
{
    // Two-site evolutions
    DmrgParameters p;
    p.set("L", 4);
    p.set("width_FermiHubbard", 2);
    p.set("height_FermiHubbard", 2);
    p.set("nsweeps", 10);
    p.set("max_bond_dimension", 10);
    p.set("u1_total_charge1", 2);
    p.set("u1_total_charge2", 1);
    p.set("symmetry", "2u1");
    p.set("MODEL", "fermi_hubbard_real");
    p.set("site_types", "0,0,0,0");
    p.set("U_FermiHubbard", 4.);
    p.set("propagator_maxiter", 10);
    p.set("imaginary_time", "yes");
    p.set("TD_backpropagation", "no");
    p.set("simulation_type", "TD");
    p.set("COMPLEX", 1);
    p.set("time_units", "fs");
    p.set("time_step", 10.);
    // iTD-DMRG
    maquis::DMRGInterface<std::complex<double>> interfaceTD(p);
    interfaceTD.evolve();
    auto energyTD = std::real(interfaceTD.energy());
    // tcDMRG
    p.set("transcorrelated_hamiltonian", "yes");
    p.set("J_Transcorrelated", 0.);
    maquis::DMRGInterface<std::complex<double>> interfaceTC(p);
    interfaceTC.evolve();
    auto energyTC = std::real(interfaceTC.energy());
    BOOST_CHECK_CLOSE(energyTD, energyTC, 1.0E-8);
}

///** @brief iTD-DMRG calculation on the 2x2 momentum-space Fermi-Hubbard model */
//BOOST_AUTO_TEST_CASE( TestImaginaryTimeFermiHubbard_2x2_Conventional_MomentumSpace )
//{
//    // Two-site evolutions
//    DmrgParameters p;
//    p.set("L", 4);
//    p.set("width_FermiHubbard", 2);
//    p.set("height_FermiHubbard", 2);
//    p.set("nsweeps", 10);
//    p.set("max_bond_dimension", 10);
//    p.set("u1_total_charge1", 2);
//    p.set("u1_total_charge2", 1);
//    p.set("symmetry", "2u1");
//    p.set("MODEL", "fermi_hubbard_momentum");
//    p.set("site_types", "0,0,0,0");
//    p.set("U_FermiHubbard", 4.);
//    p.set("propagator_maxiter", 10);
//    p.set("imaginary_time", "yes");
//    p.set("TD_backpropagation", "no");
//    p.set("simulation_type", "TD");
//    p.set("COMPLEX", 1);
//    p.set("time_units", "fs");
//    p.set("time_step", 10.);
//    // TD
//    maquis::DMRGInterface<std::complex<double>> interfaceTD(p);
//    interfaceTD.evolve();
//    auto energyTD = std::real(interfaceTD.energy());
//    // TI
//    maquis::DMRGInterface<std::complex<double>> interfaceTI(p);
//    interfaceTI.optimize();
//    auto energyTI = std::real(interfaceTI.energy());
//    BOOST_CHECK_CLOSE(energyTD, energyTI, 1.0E-8);
//}

#endif // HAVE_TwoU1