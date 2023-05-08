/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
