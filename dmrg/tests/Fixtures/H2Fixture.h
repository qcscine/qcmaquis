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
        // Real-valued integrals
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
    }
    // Class members
    RealIntegralMapType integralH2;
    DmrgParameters parametersH2;
};

#endif