/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef TEST_LiH_FIXTURE_H
#define TEST_LiH_FIXTURE_H

#include "maquis_dmrg.h"
#include "dmrg/block_matrix/symmetry.h"

/**
 * @brief Fixture class for LiH.
 * The integrals are generated for LiH with d=1 Angstrom, singlet, cc-pVDZ basis
 * set, CAS(2,4). Note that the core orbitals of Li are excluded. The data were
 * generated with MOLCAS by Leon Freitag.
 */
struct LiHFixture {
  // Types definition
  using RealIntegralMapType = typename maquis::integral_map<double>;

  /** @brief Constructor for the fixture class */
  LiHFixture() {
    // Integrals for LiH
    integralsLiH = RealIntegralMapType{{{1, 1, 1, 1}, 0.597263715971},
                                       {{1, 1, 2, 1}, 0.106899493032E-01},
                                       {{2, 1, 2, 1}, 0.215106203079E-02},
                                       {{2, 2, 2, 1}, 0.310109624236E-02},
                                       {{1, 1, 2, 2}, 0.205585394307},
                                       {{2, 2, 2, 2}, 0.266391740477},
                                       {{1, 1, 3, 1}, -0.294917730728E-08},
                                       {{2, 1, 3, 1}, 0.185509604436E-07},
                                       {{2, 2, 3, 1}, -0.437276233593E-07},
                                       {{3, 1, 3, 1}, 0.983873746561E-01},
                                       {{3, 2, 3, 1}, -0.550389961495E-02},
                                       {{3, 3, 3, 1}, -0.148187533022E-07},
                                       {{1, 1, 3, 2}, 0.656187599279E-07},
                                       {{2, 1, 3, 2}, -0.686302758417E-08},
                                       {{2, 2, 3, 2}, -0.185120456521E-07},
                                       {{3, 2, 3, 2}, 0.329840207584E-02},
                                       {{3, 3, 3, 2}, 0.705823689976E-07},
                                       {{1, 1, 3, 3}, 0.505228366944},
                                       {{2, 1, 3, 3}, 0.402391626267E-02},
                                       {{2, 2, 3, 3}, 0.207289323817},
                                       {{3, 3, 3, 3}, 0.485199794875},
                                       {{1, 1, 4, 1}, -0.179578099756},
                                       {{2, 1, 4, 1}, -0.968917151241E-02},
                                       {{2, 2, 4, 1}, -0.985834429751E-02},
                                       {{3, 1, 4, 1}, 0.166302185468E-07},
                                       {{3, 2, 4, 1}, -0.125974448467E-07},
                                       {{3, 3, 4, 1}, -0.113847869990},
                                       {{4, 1, 4, 1}, 0.118352897835},
                                       {{4, 2, 4, 1}, 0.102654021605E-01},
                                       {{4, 3, 4, 1}, -0.130590354090E-07},
                                       {{4, 4, 4, 1}, -0.121351408757},
                                       {{1, 1, 4, 2}, -0.350908680238E-01},
                                       {{2, 1, 4, 2}, 0.232966449115E-02},
                                       {{2, 2, 4, 2}, 0.137817149158E-01},
                                       {{3, 1, 4, 2}, 0.112151199425E-07},
                                       {{3, 2, 4, 2}, -0.240005477894E-07},
                                       {{3, 3, 4, 2}, -0.312450484337E-01},
                                       {{4, 2, 4, 2}, 0.123070217302E-01},
                                       {{4, 3, 4, 2}, -0.474523186140E-09},
                                       {{4, 4, 4, 2}, -0.249154148625E-01},
                                       {{1, 1, 4, 3}, 0.405834174486E-07},
                                       {{2, 1, 4, 3}, 0.416033554153E-08},
                                       {{2, 2, 4, 3}, -0.294525852258E-07},
                                       {{3, 1, 4, 3}, -0.444089181829E-02},
                                       {{3, 2, 4, 3}, -0.612942989364E-02},
                                       {{3, 3, 4, 3}, 0.350197853722E-07},
                                       {{4, 3, 4, 3}, 0.251137992170E-01},
                                       {{4, 4, 4, 3}, 0.323896546368E-07},
                                       {{1, 1, 4, 4}, 0.450421787348},
                                       {{2, 1, 4, 4}, 0.671517333359E-02},
                                       {{2, 2, 4, 4}, 0.195606443342},
                                       {{3, 1, 4, 4}, 0.185362625807E-08},
                                       {{3, 2, 4, 4}, 0.419234994857E-07},
                                       {{3, 3, 4, 4}, 0.382638069339},
                                       {{4, 4, 4, 4}, 0.370380122890},
                                       {{1, 1, 0, 0}, -0.876082926130},
                                       {{2, 1, 0, 0}, 0.928246209082E-02},
                                       {{2, 2, 0, 0}, -0.383002645306},
                                       {{3, 1, 0, 0}, 0.110478689971E-07},
                                       {{3, 2, 0, 0}, -0.838075467983E-07},
                                       {{3, 3, 0, 0}, -0.192723200850},
                                       {{4, 1, 0, 0}, 0.118275609870},
                                       {{4, 2, 0, 0}, 0.539573313879E-01},
                                       {{4, 3, 0, 0}, -0.670886115563E-07},
                                       {{4, 4, 0, 0}, -0.240135259399},
                                       {{0, 0, 0, 0}, -6.71049529388}};
    // Integrals for LiH in a STO-3G basis
    integralFileLiH_STO3GBasis.open("IntegralFile_LiH_STO3G");
    integralFileLiH_STO3GBasis << " &FCI NORB=   6,NELEC= 4,MS2=0,"
                               << std::endl;
    integralFileLiH_STO3GBasis << "  ORBSYM=1,1,1,1,1,1," << std::endl;
    integralFileLiH_STO3GBasis << "  ISYM=1," << std::endl;
    integralFileLiH_STO3GBasis << " &END" << std::endl;
    integralFileLiH_STO3GBasis
        << "1.659942303184244e+00    1         1         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.029639073398379e-01   2         1         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.049756963300113e-02    2         1         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.703227182805326e-01    2         2         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.198737934388531e-04    2         2         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "4.009795434927432e-01    2         2         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.428647316480032e-01   3         1         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.215213605109400e-02    3         1         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-7.382936572408974e-03   3         1         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.129253185855753e-02    3         1         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "6.568131760050644e-02    3         2         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-2.722016548771394e-03   3         2         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-8.953335469205959e-02   3         2         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.166943449626271e-03   3         2         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "6.103026957686852e-02    3         2         3         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.671951405644246e-01    3         3         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-6.997888264957511e-03   3         3         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.273700234025886e-01    3         3         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-9.497664232293655e-04   3         3         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.465371417146204e-02    3         3         3         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.960112231470189e-01    3         3         3         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "9.781506950959921e-03    4         1         4         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "7.759006601821306e-03    4         2         4         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.183458421130400e-02    4         2         4         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.050556899533380e-02    4         3         4         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.424221633989375e-02    4         3         4         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "4.050288847291043e-02    4         3         4         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.963524769987959e-01    4         4         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-3.577147849748006e-03   4         4         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.155942352630202e-01    4         4         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-5.030535337357253e-03   4         4         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.615973893995716e-02    4         4         3         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.663974497337714e-01    4         4         3         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.129455111594092e-01    4         4         4         4         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "9.781506950959921e-03    5         1         5         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "7.759006601821306e-03    5         2         5         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.183458421130400e-02    5         2         5         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.050556899533380e-02    5         3         5         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.424221633989375e-02    5         3         5         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "4.050288847291043e-02    5         3         5         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.686913951369104e-02    5         4         5         4         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.963524769987959e-01    5         5         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-3.577147849748006e-03   5         5         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.155942352630202e-01    5         5         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-5.030535337357253e-03   5         5         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.615973893995716e-02    5         5         3         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.663974497337714e-01    5         5         3         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.792072321320272e-01    5         5         4         4         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.129455111594092e-01    5         5         5         5         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "5.021536562889917e-02    6         1         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-7.107540925710492e-03   6         1         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-5.902086725461253e-03   6         1         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-2.562737813835795e-03   6         1         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.249990153850549e-03    6         1         3         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "9.955159245318470e-03    6         1         3         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.327853106968450e-03    6         1         4         4         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.327853106968450e-03    6         1         5         5         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "9.260399182385273e-03    6         1         6         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-9.128541548956752e-02   6         2         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.535223351372107e-04    6         2         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "9.111393914190066e-02    6         2         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "5.177792680017361e-03    6         2         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-7.339950754819091e-02   6         2         3         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.399648874423638e-03    6         2         3         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-4.940584542147215e-02   6         2         4         4         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-4.940584542147215e-02   6         2         5         5         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.618750434179279e-03    6         2         6         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.215936835174049e-01    6         2         6         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "4.331064128286324e-02    6         3         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-2.278155107173152e-03   6         3         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-8.145293602036821e-02   6         3         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.668633235939294e-03    6         3         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "4.998493455473627e-02    6         3         3         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.122485961024775e-02    6         3         3         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.188298536918154e-02    6         3         4         4         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.188298536918154e-02    6         3         5         5         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "6.370511344668381e-03    6         3         6         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-5.185367132945327e-02   6         3         6         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "5.824935217396741e-02    6         3         6         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-4.095031188959358e-03   6         4         4         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.455528664518197e-02   6         4         4         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-6.840852841711873e-03   6         4         4         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.658528723767585e-02    6         4         6         4         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-4.095031188959358e-03   6         5         5         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.455528664518197e-02   6         5         5         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-6.840852841711873e-03   6         5         5         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.658528723767585e-02    6         5         6         5         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.423343879861308e-01    6         6         1         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-9.209922072783469e-04   6         6         2         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.481692451442813e-01    6         6         2         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-8.161717960324199e-03   6         6         3         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-4.699417976310499e-02   6         6         3         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.521057175798820e-01    6         6         3         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.496314991699705e-01    6         6         4         4         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "2.496314991699705e-01    6         6         5         5         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-5.049014224854249e-03   6         6         6         1         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.555853752941223e-02    6         6         6         2         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-4.149504401053612e-02   6         6         6         3         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.377252763254533e-01    6         6         6         6         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-4.573998049531848e+00   1         1         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.028440336841015e-01    2         1         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.106614296303322e+00   2         2         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "1.549085896837615e-01    3         1         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-2.967715708861850e-02   3         2         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.049578154428433e+00   3         3         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.041179356238942e+00   4         4         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.041179356238942e+00   5         5         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-3.815766726733125e-02   6         1         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "8.434934497920570e-02    6         2         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "3.223434847057938e-04    6         3         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "-1.015815175560325e+00   6         6         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis
        << "5.291772109200000e-01    0         0         0         0         "
        << std::endl;
    integralFileLiH_STO3GBasis.close();

    // LiH parameters
    // NOTE: Is this really LiH, nelec=2 looks suspicious
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

    // LiH parameters for the STO3G basis
    parametersLiH_STO3G.set("integral_file", "IntegralFile_LiH_STO3G");
    parametersLiH_STO3G.set("L", 6);
    parametersLiH_STO3G.set("irrep", 0);
    parametersLiH_STO3G.set("nsweeps", 2);
    parametersLiH_STO3G.set("max_bond_dimension", 100);
    parametersLiH_STO3G.set("nelec", 4);
    parametersLiH_STO3G.set("spin", 0);
    parametersLiH_STO3G.set("u1_total_charge1", 2);
    parametersLiH_STO3G.set("u1_total_charge2", 2);
    parametersLiH_STO3G.set("MODEL", "quantum_chemistry");

    // CISD wave function for LiH in a minimal basis (obtained from PySCF)
    // 4,4,1,1,1,1       -0.8487587500849899
    // 4,1,4,1,1,1       0.29026577890912375
    // 4,1,1,4,1,1       0.010372012270461816
    // 4,1,1,1,4,1       0.010372012270461816
    // 4,1,1,1,1,4       0.176659694841801
    // 1,4,4,1,1,1       0.0019936569973847476
    // 1,4,1,4,1,1       0.0015602133654161537
    // 1,4,1,1,4,1       0.0015602133654161537
    // 1,4,1,1,1,4       0.0011716368962185644
    //
    // 4,3,2,1,1,1       -0.15525689943121151
    // 4,2,3,1,1,1       0.15525689943121151
    // 4,3,1,1,1,2       -0.04267287376249823
    // 4,2,1,1,1,3       0.04267287376249823
    // 4,1,3,1,1,2       -0.23662594098756534
    // 4,1,2,1,1,3       0.23662594098756534
    // 3,4,2,1,1,1       -0.0001190762369702869
    // 2,4,3,1,1,1       0.0001190762369702869
    // 3,4,1,1,1,2       0.0011686475234696796
    // 2,4,1,1,1,3       -0.0011686475234696796
    // 3,2,1,1,1,4       -0.002061401722547903
    // 2,3,1,1,1,4       0.0020614017225479026
    // 3,2,4,1,1,1       0.0006293871978559377
    // 2,3,4,1,1,1       -0.0006293871978559378
    // 1,4,3,1,1,2       0.0002559946198255345
    // 1,4,2,1,1,3       -0.000255994619825532
    // 3,2,1,4,1,1       -0.0034073198652948395
    // 2,3,1,4,1,1       0.00340731986529484
    // 3,2,1,1,4,1       -0.0034073198652948395
    // 2,3,1,1,4,1       0.0034073198652948
    //
    // 3,2,2,1,1,3       0.001991837144718305
    // 2,3,2,1,1,3       -0.0029186405027821685
    // 2,2,3,1,1,3       0.0009268033580638636
    // 3,3,2,1,1,2       0.0009268033580638636
    // 3,2,3,1,1,2       -0.002918640502782168
    // 2,3,3,1,1,2       0.001991837144718305

    cisdWavefunctionLiHSto3g =
        std::string(
            "4,4,1,1,1,1|4,1,4,1,1,1|4,1,1,4,1,1|4,1,1,1,4,1|4,1,1,1,1,4|"
        ) +
        std::string("1,4,4,1,1,1|1,4,1,4,1,1|1,4,1,1,4,1|1,4,1,1,1,4|") +
        std::string("4,3,2,1,1,1|4,2,3,1,1,1|4,3,1,1,1,2|4,2,1,1,1,3|") +
        std::string("4,1,3,1,1,2|4,1,2,1,1,3|3,4,2,1,1,1|2,4,3,1,1,1|") +
        std::string("3,4,1,1,1,2|2,4,1,1,1,3|3,2,1,1,1,4|2,3,1,1,1,4|") +
        std::string("3,2,4,1,1,1|2,3,4,1,1,1|1,4,3,1,1,2|1,4,2,1,1,3|") +
        std::string("3,2,1,4,1,1|2,3,1,4,1,1|3,2,1,1,4,1|2,3,1,1,4,1");
    //+
    // std::string("3,2,2,1,1,3|2,3,2,1,1,3|2,2,3,1,1,3|3,3,2,1,1,2|3,2,3,1,1,2|2,3,3,1,1,2");

    // cisdWavefunctionLiHSto3g =
    // std::string("4,4,1,1,1,1|4,3,2,1,1,1|3,4,2,1,1,1|4,3,1,1,1,2|3,4,1,1,1,2|")
    //                          +
    //                          std::string("3,3,2,1,1,2|4,2,3,1,1,1|4,1,4,1,1,1|3,2,4,1,1,1|4,1,3,1,1,2|")
    //                          +
    //                          std::string("3,2,3,1,1,2|2,4,3,1,1,1|2,3,4,1,1,1|1,4,4,1,1,1|2,3,3,1,1,2|")
    //                          +
    //                          std::string("1,4,3,1,1,2|4,1,1,4,1,1|3,2,1,4,1,1|2,3,1,4,1,1|1,4,1,4,1,1|")
    //                          +
    //                          std::string("4,1,1,1,4,1|3,2,1,1,4,1|2,3,1,1,4,1|1,4,1,1,4,1|4,2,1,1,1,3|")
    //                          +
    //                          std::string("4,1,2,1,1,3|3,2,2,1,1,3|4,1,1,1,1,4|3,2,1,1,1,4|2,4,1,1,1,3|")
    //                          +
    //                          std::string("2,3,2,1,1,3|1,4,2,1,1,3|2,3,1,1,1,4|1,4,1,1,1,4|2,2,3,1,1,3");

    cisdCoefficientsLiHSto3g =
        std::string(
            "-0.8487587500849899,0.29026577890912375,0.01037201227046181,0."
            "010372012270461816,0.176659694841801,"
        ) +
        std::string(
            "0.0019936569973847476,0.0015602133654161537,0.0015602133654161537,"
            "0.0011716368962185644,"
        ) +
        std::string(
            "-0.15525689943121151,0.15525689943121151,-0.04267287376249823,0."
            "04267287376249823,"
        ) +
        std::string(
            "-0.23662594098756534,0.23662594098756534,-0.0001190762369702869,0."
            "0001190762369702869,"
        ) +
        std::string(
            "0.0011686475234696796,-0.0011686475234696796,-0."
            "002061401722547903,0.0020614017225479026,"
        ) +
        std::string(
            "0.0006293871978559377,-0.0006293871978559378,0."
            "0002559946198255345,-0.000255994619825532,"
        ) +
        std::string(
            "-0.0034073198652948395,0.00340731986529484,-0.0034073198652948395,"
            "0.00340731986529484"
        );
    // +
    // std::string("0.001991837144718305,-0.0029186405027821685,0.0009268033580638636,0.0009268033580638636,-0.002918640502782168,0.001991837144718305");

    // cisdCoefficientsLiHSto3g =
    // std::string("-0.8487587500849899,-0.15525689943121151,-0.0001190762369702869,")
    //                          +
    //                          std::string("-0.04267287376249823,0.0011686475234696796,0.0009268033580638636,")
    //                          +
    //                          std::string("0.15525689943121151,0.29026577890912375,0.0006293871978559377,")
    //                          +
    //                          std::string("-0.23662594098756534,-0.002918640502782168,0.0001190762369702869,")
    //                          +
    //                          std::string("-0.0006293871978559378,0.0019936569973847476,0.001991837144718305,")
    //                          +
    //                          std::string("0.0002559946198255345,0.010372012270461816,-0.0034073198652948395,")
    //                          +
    //                          std::string("0.00340731986529484,0.0015602133654161537,0.010372012270461816,")
    //                          +
    //                          std::string("-0.0034073198652948395,0.00340731986529484,0.0015602133654161537,")
    //                          +
    //                          std::string("0.04267287376249823,0.23662594098756534,0.001991837144718305,")
    //                          +
    //                          std::string("0.176659694841801,-0.002061401722547903,-0.0011686475234696796,")
    //                          +
    //                          std::string("-0.0029186405027821685,-0.000255994619825532,0.0020614017225479026,")
    //                          +
    //                          std::string("0.0011716368962185644,0.0009268033580638636");

    cisdWavefunctionLiHSto3gNR =
        std::string(
            "4,4,1,1,1,1|4,1,4,1,1,1|4,1,1,4,1,1|4,1,1,1,4,1|4,1,1,1,1,4|"
        ) +
        std::string("1,4,4,1,1,1|1,4,1,4,1,1|1,4,1,1,4,1|1,4,1,1,1,4|") +
        std::string("4,3,2,1,1,1|4,3,1,1,1,2|") +
        std::string("4,1,3,1,1,2|3,4,2,1,1,1|") +
        std::string("3,4,1,1,1,2|3,2,1,1,1,4|") +
        std::string("3,2,4,1,1,1|1,4,3,1,1,2|") +
        std::string("3,2,1,4,1,1|3,2,1,1,4,1");
    //+ std::string("3,3,2,1,1,2");

    // cisdWavefunctionLiHSto3gNR =
    // std::string("4,4,1,1,1,1|4,3,2,1,1,1|3,4,2,1,1,1|4,3,1,1,1,2|3,4,1,1,1,2|")
    //                            +
    //                            std::string("4,1,4,1,1,1|3,2,4,1,1,1|4,1,3,1,1,2|1,4,4,1,1,1|1,4,3,1,1,2|")
    //                            +
    //                            std::string("4,1,1,4,1,1|3,2,1,4,1,1|1,4,1,4,1,1|4,1,1,1,4,1|3,2,1,1,4,1|")
    //                            +
    //                            std::string("1,4,1,1,4,1|4,1,1,1,1,4|3,2,1,1,1,4|1,4,1,1,1,4|")
    //                            +
    //                            std::string("3,3,2,1,1,2|3,2,3,1,1,2|2,3,3,1,1,2|3,2,2,1,1,3|2,3,2,1,1,3|2,2,3,1,1,3");

    cisdCoefficientsLiHSto3gNR =
        std::string(
            "-0.8487587500849899,0.29026577890912375,0.01037201227046181,0."
            "010372012270461816,0.176659694841801,"
        ) +
        std::string(
            "0.0019936569973847476,0.0015602133654161537,0.0015602133654161537,"
            "0.0011716368962185644,"
        ) +
        std::string("-0.219566412827615,-0.06034855682036001,") +
        std::string("-0.3346396149539106,-0.00016839922927973228,") +
        std::string("0.0016527171773245508,-0.0029152622735265045,") +
        std::string("0.0008900879111918658,0.0003620310632518153,") +
        std::string("-0.00481867796484323,-0.00481867796484323");
    //+ std::string("0.002270195319154438");

    // cisdCoefficientsLiHSto3gNR =
    // std::string("-0.8487587500849899,-0.219566412827615,-0.00016839922927973228,-0.06034855682036001,0.0016527171773245508,")
    //                            +
    //                            std::string("0.29026577890912375,0.0008900879111918658,-0.3346396149539106,0.0019936569973847476,0.0003620310632518153,")
    //                            +
    //                            std::string("0.010372012270461816,-0.00481867796484323,0.0015602133654161537,0.010372012270461816,-0.00481867796484323,")
    //                            +
    //                            std::string("0.0022064749015672943,0.176659694841801,-0.0029152622735265045,0.0011716368962185644,")
    //                            +
    //                            std::string("0.002270195319154439,-0.007149179974436458,0.00487898465528202,0.00487898465528202,-0.007149179974436458,0.002270195319154439");
  }

  /** @brief Removes the FCIDUMP file */
  ~LiHFixture() { std::remove("IntegralFile_LiH_STO3G"); }

  // Class members
  RealIntegralMapType integralsLiH;
  DmrgParameters parametersLiH, parametersLiH_STO3G;
  double referenceEnergy =
      -7.90435750473166;  // This reference is taken from test2.cpp
  std::ofstream integralFileLiH_STO3GBasis;
  std::string cisdWavefunctionLiHSto3g, cisdCoefficientsLiHSto3g,
      cisdWavefunctionLiHSto3gNR, cisdCoefficientsLiHSto3gNR;
  double referenceSto3gEnergy = -7.798753320072483,
         referenceSto3gHFEnergy = -7.71082990021723;
};

#endif
