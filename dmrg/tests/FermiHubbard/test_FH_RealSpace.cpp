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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include "utils/io.hpp"
#include <iostream>
#include "maquis_dmrg.h" // Needed for the interface

/**
 * @brief Test on the energy of a 2x2 real-space Fermi-Hubbard model
 * Reference energy taken from "Symmetry in auxiliary-field quantum Monte 
 * Carlo calculations", PRB, 2013
 */
BOOST_AUTO_TEST_CASE( Test_FermiHubbardRealSpace_2x2_2Alpha1Beta )
{
#ifdef HAVE_TwoU1
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
    maquis::DMRGInterface<double> interface(p);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), -1.6046*4, 1.0e-2);
#endif // HAVE_TwoU1
}

#ifdef HAVE_TwoU1

/** @brief Test on the energy of a 3x3 real-space Fermi-Hubbard model */
BOOST_AUTO_TEST_CASE( Test_FermiHubbardRealSpace_3x3_4Alpha4Beta )
{
    DmrgParameters p;
    p.set("L", 9);
    p.set("width_FermiHubbard", 3);
    p.set("height_FermiHubbard", 3);
    p.set("nsweeps", 10);
    p.set("max_bond_dimension", 300);
    p.set("u1_total_charge1", 4);
    p.set("u1_total_charge2", 4);
    p.set("symmetry", "2u1");
    p.set("MODEL", "fermi_hubbard_real");
    p.set("site_types", "0,0,0,0,0,0,0,0,0");
    p.set("U_FermiHubbard", 8.);
    maquis::DMRGInterface<double> interface(p);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), -0.8094*9, 1.0e-2);
}

#endif // HAVE_TwoU1