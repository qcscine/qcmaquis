/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef TEST_LiH_FIXTURE_H
#define TEST_LiH_FIXTURE_H

#include "maquis_dmrg.h"
#include "dmrg/block_matrix/symmetry.h"

/**
 * @brief Fixture class for LiH.
 * The integrals are generated for LiH with d=1 Angstrom, singlet, cc-pVDZ basis set, CAS(2,4).
 * Note that the core orbitals of Li are excluded.
 * The data were generated with MOLCAS by Leon Freitag.
 */
struct LiHFixture
{
// Types definition
  using RealIntegralMapType = typename maquis::integral_map<double>;

  /** @brief Constructor for the fixture class */
  LiHFixture() {
    // Real-valued integrals
    integralsLiH = RealIntegralMapType {
                                        { { 1, 1, 1, 1 },   0.597263715971     },
                                        { { 1, 1, 2, 1 },   0.106899493032E-01 },
                                        { { 2, 1, 2, 1 },   0.215106203079E-02 },
                                        { { 2, 2, 2, 1 },   0.310109624236E-02 },
                                        { { 1, 1, 2, 2 },   0.205585394307     },
                                        { { 2, 2, 2, 2 },   0.266391740477     },
                                        { { 1, 1, 3, 1 },  -0.294917730728E-08 },
                                        { { 2, 1, 3, 1 },   0.185509604436E-07 },
                                        { { 2, 2, 3, 1 },  -0.437276233593E-07 },
                                        { { 3, 1, 3, 1 },   0.983873746561E-01 },
                                        { { 3, 2, 3, 1 },  -0.550389961495E-02 },
                                        { { 3, 3, 3, 1 },  -0.148187533022E-07 },
                                        { { 1, 1, 3, 2 },   0.656187599279E-07 },
                                        { { 2, 1, 3, 2 },  -0.686302758417E-08 },
                                        { { 2, 2, 3, 2 },  -0.185120456521E-07 },
                                        { { 3, 2, 3, 2 },   0.329840207584E-02 },
                                        { { 3, 3, 3, 2 },   0.705823689976E-07 },
                                        { { 1, 1, 3, 3 },   0.505228366944     },
                                        { { 2, 1, 3, 3 },   0.402391626267E-02 },
                                        { { 2, 2, 3, 3 },   0.207289323817     },
                                        { { 3, 3, 3, 3 },   0.485199794875     },
                                        { { 1, 1, 4, 1 },  -0.179578099756     },
                                        { { 2, 1, 4, 1 },  -0.968917151241E-02 },
                                        { { 2, 2, 4, 1 },  -0.985834429751E-02 },
                                        { { 3, 1, 4, 1 },   0.166302185468E-07 },
                                        { { 3, 2, 4, 1 },  -0.125974448467E-07 },
                                        { { 3, 3, 4, 1 },  -0.113847869990     },
                                        { { 4, 1, 4, 1 },   0.118352897835     },
                                        { { 4, 2, 4, 1 },   0.102654021605E-01 },
                                        { { 4, 3, 4, 1 },  -0.130590354090E-07 },
                                        { { 4, 4, 4, 1 },  -0.121351408757     },
                                        { { 1, 1, 4, 2 },  -0.350908680238E-01 },
                                        { { 2, 1, 4, 2 },   0.232966449115E-02 },
                                        { { 2, 2, 4, 2 },   0.137817149158E-01 },
                                        { { 3, 1, 4, 2 },   0.112151199425E-07 },
                                        { { 3, 2, 4, 2 },  -0.240005477894E-07 },
                                        { { 3, 3, 4, 2 },  -0.312450484337E-01 },
                                        { { 4, 2, 4, 2 },   0.123070217302E-01 },
                                        { { 4, 3, 4, 2 },  -0.474523186140E-09 },
                                        { { 4, 4, 4, 2 },  -0.249154148625E-01 },
                                        { { 1, 1, 4, 3 },   0.405834174486E-07 },
                                        { { 2, 1, 4, 3 },   0.416033554153E-08 },
                                        { { 2, 2, 4, 3 },  -0.294525852258E-07 },
                                        { { 3, 1, 4, 3 },  -0.444089181829E-02 },
                                        { { 3, 2, 4, 3 },  -0.612942989364E-02 },
                                        { { 3, 3, 4, 3 },   0.350197853722E-07 },
                                        { { 4, 3, 4, 3 },   0.251137992170E-01 },
                                        { { 4, 4, 4, 3 },   0.323896546368E-07 },
                                        { { 1, 1, 4, 4 },   0.450421787348     },
                                        { { 2, 1, 4, 4 },   0.671517333359E-02 },
                                        { { 2, 2, 4, 4 },   0.195606443342     },
                                        { { 3, 1, 4, 4 },   0.185362625807E-08 },
                                        { { 3, 2, 4, 4 },   0.419234994857E-07 },
                                        { { 3, 3, 4, 4 },   0.382638069339     },
                                        { { 4, 4, 4, 4 },   0.370380122890     },
                                        { { 1, 1, 0, 0 },  -0.876082926130     },
                                        { { 2, 1, 0, 0 },   0.928246209082E-02 },
                                        { { 2, 2, 0, 0 },  -0.383002645306     },
                                        { { 3, 1, 0, 0 },   0.110478689971E-07 },
                                        { { 3, 2, 0, 0 },  -0.838075467983E-07 },
                                        { { 3, 3, 0, 0 },  -0.192723200850     },
                                        { { 4, 1, 0, 0 },   0.118275609870     },
                                        { { 4, 2, 0, 0 },   0.539573313879E-01 },
                                        { { 4, 3, 0, 0 },  -0.670886115563E-07 },
                                        { { 4, 4, 0, 0 },  -0.240135259399     },
                                        { { 0, 0, 0, 0 },   -6.71049529388     }
                                      };
    // LiH parameters
    parametersLiH.set("integrals_binary", maquis::serialize(integralsLiH));
    parametersLiH.set("site_types", "0,0,0,0");
    parametersLiH.set("L", 4);
    parametersLiH.set("irrep", 0);
    parametersLiH.set("nsweeps", 2);
    parametersLiH.set("max_bond_dimension", 100);
    parametersLiH.set("nelec", 2);
    parametersLiH.set("spin", 0);
    parametersLiH.set("u1_total_charge1", 1);
    parametersLiH.set("u1_total_charge2", 1);
    parametersLiH.set("MODEL", "quantum_chemistry");
  }
  // Class members
  RealIntegralMapType integralsLiH;
  DmrgParameters parametersLiH;
  double referenceEnergy = -7.90435750473166; // This reference is taken from test2.cpp
};

#endif
