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

#ifndef TEST_H2_FIXTURE_H
#define TEST_H2_FIXTURE_H

#include "maquis_dmrg.h"
#include "dmrg/block_matrix/symmetry.h"

/**
 * @brief Fixture class for H2.
 */
struct H2Fixture
{
    // Types definition
    using RealIntegralMapType = typename maquis::integral_map<double>;

    /** @brief Constructor for the fixture class */
    H2Fixture() {
        // H2 integrals in the conventional format
        integralH2 = RealIntegralMapType {
                                            { { 1, 1, 1, 1 },   0.354237848011       },
                                            { { 1, 1, 2, 1 },  -0.821703816101E-13   },
                                            { { 2, 1, 2, 1 },  0.185125251547        },
                                            { { 2, 2, 2, 1 },  0.782984788117E-13    },
                                            { { 1, 1, 2, 2 },  0.361001163519        },
                                            { { 2, 2, 2, 2 },  0.371320200119        },
                                            { { 1, 1, 0, 0 }, -0.678487901790        },
                                            { { 2, 1, 0, 0 }, -0.539801158857E-14    },
                                            { { 2, 2, 0, 0 }, -0.653221638776        },
                                            { { 0, 0, 0, 0 },  0.176392403557        }
                                          };

        //  H2 integral file in conventional format
        integralFileH2ConventionalFormat.open("IntegralFile_H2_ConventionalFormat");
        integralFileH2ConventionalFormat << " &FCI NORB=2,NELEC=2,MS2=0, " << std::endl;
        integralFileH2ConventionalFormat << " ORBSYM=1,1," << std::endl;
        integralFileH2ConventionalFormat << " ISYM=1, " << std::endl;
        integralFileH2ConventionalFormat << " &END " << std::endl;
        integralFileH2ConventionalFormat << " 0.354237848011E+00   1  1  1  1" << std::endl;
        integralFileH2ConventionalFormat << " 0.185125251547E+00   2  1  2  1" << std::endl;
        integralFileH2ConventionalFormat << " 0.361001163519E+00   1  1  2  2" << std::endl;
        integralFileH2ConventionalFormat << " 0.371320200119E+00   2  2  2  2" << std::endl;
        integralFileH2ConventionalFormat << "-0.678487901790E+00   1  1  0  0" << std::endl;
        integralFileH2ConventionalFormat << "-0.653221638776E+00   2  2  0  0" << std::endl;
        integralFileH2ConventionalFormat << " 0.176392403557E+00   0  0  0  0" << std::endl;
        integralFileH2ConventionalFormat.close();

        // H2 integral file in quantum format.
        integralFileH2QuantumFormat.open("IntegralFile_H2_QuantumFormat");
        integralFileH2QuantumFormat << " &FCI NORB=2,NELEC=2,MS2=0, " << std::endl;
        integralFileH2QuantumFormat << " ORBSYM=1,1," << std::endl;
        integralFileH2QuantumFormat << " ISYM=1, " << std::endl;
        integralFileH2QuantumFormat << " &END " << std::endl;
        integralFileH2QuantumFormat << "-9.481694515690000e-01    1   1    0    0    0    0" << std::endl;
        integralFileH2QuantumFormat << "-9.314443646090000e-01    2   2    0    0    0    0" << std::endl;
        integralFileH2QuantumFormat << "+3.542378480110000e-01    1   1    1    1    0    0" << std::endl;
        integralFileH2QuantumFormat << "+1.851252515470000e-01    1   2    1    2    0    0" << std::endl;
        integralFileH2QuantumFormat << "+1.851252515470000e-01    2   1    1    2    0    0" << std::endl;
        integralFileH2QuantumFormat << "+1.851252515470000e-01    2   1    2    1    0    0" << std::endl;
        integralFileH2QuantumFormat << "+3.610011635190000e-01    2   2    1    1    0    0" << std::endl;
        integralFileH2QuantumFormat << "+3.713202001190000e-01    2   2    2    2    0    0" << std::endl;
        integralFileH2QuantumFormat << "+0.176392403557000E+00    0   0    0    0    0    0" << std::endl;
        integralFileH2QuantumFormat.close();

        // == Conventional calculation ==
        // Generic parameters
        parametersH2.set("integrals_binary", maquis::serialize(integralH2));
        parametersH2.set("site_types", "0,0");
        parametersH2.set("L", 2);
        parametersH2.set("irrep", 0);
        parametersH2.set("nsweeps",2);
        parametersH2.set("max_bond_dimension",100);
        // for SU2U1
        parametersH2.set("nelec", 2);
        parametersH2.set("spin", 0);
        // for 2U1
        parametersH2.set("u1_total_charge1", 1);
        parametersH2.set("u1_total_charge2", 1);

        // == Quantum format ==
        // -- General parameters --
        parametersH2QuantumFormat.set("L", 2);
        parametersH2QuantumFormat.set("max_bond_dimension", 10);
        parametersH2QuantumFormat.set("nsweeps", 10);
        parametersH2QuantumFormat.set("symmetry", "2u1pg");
        parametersH2QuantumFormat.set("u1_total_charge1", 1);
        parametersH2QuantumFormat.set("u1_total_charge2", 1);
        parametersH2QuantumFormat.set("LATTICE", "orbitals");
        parametersH2QuantumFormat.set("CONSERVED_QUANTUMNUMBERS", "Nup,Ndown");
        parametersH2QuantumFormat.set("MODEL", "quantum_chemistry");
        parametersH2QuantumFormat.set("integral_file", "IntegralFile_H2_ConventionalFormat");
        parametersH2QuantumFormat.set("optimization", "twosite");
        parametersH2QuantumFormat.set("init_type", "default");
        parametersH2QuantumFormat.set("seed", "16071991");
        // -- Transcorrelated options --
        parametersH2QuantumFormat.set("transcorrelated_nsweeps_TI", 0);
        parametersH2QuantumFormat.set("transcorrelated_nsweeps_TC", 100);
        parametersH2QuantumFormat.set("transcorrelated_integral_file", "IntegralFile_H2_QuantumFormat");
        parametersH2QuantumFormat.set("transcorrelated_quantum_computing_format", "yes");
        // -- Time-dependent simulation --
        parametersH2QuantumFormat.set("propagator_accuracy", 1.0E-5);
        parametersH2QuantumFormat.set("propagator_maxiter", 10);
        parametersH2QuantumFormat.set("time_step", 1);
        parametersH2QuantumFormat.set("TD_backpropagation", "no");
    }

    /** @brief Class destructor */
    ~H2Fixture() {
        std::remove("IntegralFile_H2_ConventionalFormat");
        std::remove("IntegralFile_H2_QuantumFormat");
    }
    // Class members
    RealIntegralMapType integralH2;
    DmrgParameters parametersH2, parametersH2QuantumFormat;
    std::ofstream integralFileH2ConventionalFormat, integralFileH2QuantumFormat;
};

#endif
