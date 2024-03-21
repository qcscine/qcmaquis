/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#ifndef TEST_NMODE_FIXTURE_H
#define TEST_NMODE_FIXTURE_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

/**
 * @brief Fixture class for the test of the n-mode vibrational DMRG code.
 */
struct NModeFixture {
  // Types definition
  using MaquisIntegralType =
      maquis::integral_map<double, chem::Hamiltonian::VibrationalNMode>;
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
    // == PARAMETERS FOR THE ONE-MODE FAD CALCULATION ==
    // These data were generated based on the PES published by Bowman and based
    // on a DVR basis set.
    parametersFADOneBody.set("L", 39);
    parametersFADOneBody.set("nmode_num_modes", 1);
    parametersFADOneBody.set("nmode_max_coupling", 1);
    parametersFADOneBody.set("nmode_num_basis", "39");
    parametersFADOneBody.set("symmetry", "nu1");
    parametersFADOneBody.set("LATTICE", "nmode lattice");
    parametersFADOneBody.set("MODEL", "nmode");
    parametersFADOneBody.set("integral_file", "integral_file_test_OneBodyFAD");
    //
    parametersFADOneBodyPaired.set("L", 1);
    parametersFADOneBodyPaired.set("nmode_num_modes", 1);
    parametersFADOneBodyPaired.set("nmode_max_coupling", 1);
    parametersFADOneBodyPaired.set("nmode_num_basis", "39");
    parametersFADOneBodyPaired.set("symmetry", "none");
    parametersFADOneBodyPaired.set("LATTICE", "watson lattice");
    parametersFADOneBodyPaired.set("MODEL", "nmodecompactpaired");
    parametersFADOneBodyPaired.set("integral_file", "integral_file_test_OneBodyFAD");
    //
    parametersFADTwoBody.set("L", 22);
    parametersFADTwoBody.set("nmode_num_modes", 2);
    parametersFADTwoBody.set("nmode_max_coupling", 2);
    parametersFADTwoBody.set("nmode_num_basis", "11,11");
    parametersFADTwoBody.set("symmetry", "nu1");
    parametersFADTwoBody.set("LATTICE", "nmode lattice");
    parametersFADTwoBody.set("MODEL", "nmode");
    parametersFADTwoBody.set("integral_file", "integral_file_test_TwoBodyFAD");
    //
    parametersWater.set("nsweeps", 50);
    parametersWater.set("ngrowsweeps", 10);
    parametersWater.set("nmainsweeps", 20);
    parametersWater.set("max_bond_dimension", 150);
    parametersWater.set("alpha_initial", 1.0E-8);
    parametersWater.set("alpha_main", 1.0E-10);
    parametersWater.set("alpha_final", 0);
    parametersWater.set("truncation_initial", 0);
    parametersWater.set("truncation_final", 0);
    parametersWater.set("eigensolver", "IETL_JCD");
    parametersWater.set("optimization", "singlesite");
    parametersWater.set("model_library", "coded");
    parametersWater.set("lattice_library", "coded");
    parametersWater.set("integral_cutoff", 1.0E-8);
    parametersWater.set("nmode_num_modes", 3);
    parametersWater.set("init_type", "basis_state_generic");

    // == INPUT FILE CREATIONS ==
    integralFileOneBodyFAD.open("integral_file_test_OneBodyFAD");
    integralFileOneBodyFAD << "       1 0       1 0  -2.359242429009664e+03\n";
    integralFileOneBodyFAD << "       1 1       1 1  -2.358797386438359e+03\n";
    integralFileOneBodyFAD << "       1 2       1 2  -1.444541233850135e+03\n";
    integralFileOneBodyFAD << "       1 3       1 3  -1.437160009122911e+03\n";
    integralFileOneBodyFAD << "       1 4       1 4  -6.657064784181758e+03\n";
    integralFileOneBodyFAD << "       1 5       1 5  -6.128584661601011e+03\n";
    integralFileOneBodyFAD << "       1 6       1 6  -3.441604675890336e+03\n";
    integralFileOneBodyFAD << "       1 7       1 7   1.538026675462226e+03\n";
    integralFileOneBodyFAD << "       1 8       1 8   6.026905892315631e+03\n";
    integralFileOneBodyFAD << "       1 9       1 9   9.519079627438325e+03\n";
    integralFileOneBodyFAD << "      1 10      1 10   1.408367346423375e+03\n";
    integralFileOneBodyFAD << "      1 11      1 11   1.877962833530346e+03\n";
    integralFileOneBodyFAD << "      1 12      1 12   2.406532356803376e+03\n";
    integralFileOneBodyFAD << "      1 13      1 13   2.974933802551137e+03\n";
    integralFileOneBodyFAD << "      1 14      1 14   3.590607405365762e+03\n";
    integralFileOneBodyFAD << "      1 15      1 15   4.250019051366812e+03\n";
    integralFileOneBodyFAD << "      1 16      1 16   4.954168295956210e+03\n";
    integralFileOneBodyFAD << "      1 17      1 17   5.702284980641777e+03\n";
    integralFileOneBodyFAD << "      1 18      1 18   6.494407879568767e+03\n";
    integralFileOneBodyFAD << "      1 19      1 19   7.330063654159827e+03\n";
    integralFileOneBodyFAD << "      1 20      1 20   8.209556409406479e+03\n";
    integralFileOneBodyFAD << "      1 21      1 21   9.132461752263011e+03\n";
    integralFileOneBodyFAD << "      1 22      1 22   1.009901501528866e+03\n";
    integralFileOneBodyFAD << "      1 23      1 23   1.110869147300428e+03\n";
    integralFileOneBodyFAD << "      1 24      1 24   1.216206565203700e+03\n";
    integralFileOneBodyFAD << "      1 25      1 25   1.325854103259667e+03\n";
    integralFileOneBodyFAD << "      1 26      1 26   1.439864696031767e+03\n";
    integralFileOneBodyFAD << "      1 27      1 27   1.558157142548529e+03\n";
    integralFileOneBodyFAD << "      1 28      1 28   1.680827752856561e+03\n";
    integralFileOneBodyFAD << "      1 29      1 29   1.807804326901167e+03\n";
    integralFileOneBodyFAD << "      1 30      1 30   1.939221830929778e+03\n";
    integralFileOneBodyFAD << "      1 31      1 31   2.074605464868975e+03\n";
    integralFileOneBodyFAD << "      1 32      1 32   2.214735474735023e+03\n";
    integralFileOneBodyFAD << "      1 33      1 33   2.359908634757324e+03\n";
    integralFileOneBodyFAD << "      1 34      1 34   2.503011884645919e+03\n";
    integralFileOneBodyFAD << "      1 35      1 35   2.651668901080508e+03\n";
    integralFileOneBodyFAD << "      1 36      1 36   2.861283145084534e+03\n";
    integralFileOneBodyFAD << "      1 37      1 37   2.946541753544523e+03\n";
    integralFileOneBodyFAD << "      1 37      1 38   2.795969730154564e-13\n";
    integralFileOneBodyFAD << "      1 38       1 0  -4.721452208883015e-13\n";
    integralFileOneBodyFAD << "      1 38       1 1   9.167304342968547e-13\n";
    integralFileOneBodyFAD << "      1 38       1 2   3.627696363978469e-13\n";
    integralFileOneBodyFAD << "      1 38       1 3   2.368394571952736e-13\n";
    integralFileOneBodyFAD << "      1 38       1 4  -1.389061557427853e-13\n";
    integralFileOneBodyFAD << "      1 38       1 5   1.503577338928864e-13\n";
    integralFileOneBodyFAD << "      1 38       1 6  -3.152902036131743e-13\n";
    integralFileOneBodyFAD << "      1 38       1 7   2.076156246433921e-13\n";
    integralFileOneBodyFAD << "      1 38       1 8   2.296264501786515e-13\n";
    integralFileOneBodyFAD << "      1 38       1 9   8.566375343452281e-13\n";
    integralFileOneBodyFAD << "      1 38      1 10   2.391446450047095e-13\n";
    integralFileOneBodyFAD << "      1 38      1 11   2.260571271188797e-13\n";
    integralFileOneBodyFAD << "      1 38      1 12   6.781713813566390e-13\n";
    integralFileOneBodyFAD << "      1 38      1 13   2.974435883143153e-13\n";
    integralFileOneBodyFAD << "      1 38      1 14  -1.998820913472199e-13\n";
    integralFileOneBodyFAD << "      1 38      1 15   4.568733516507884e-13\n";
    integralFileOneBodyFAD << "      1 38      1 16   3.331368189120332e-13\n";
    integralFileOneBodyFAD << "      1 38      1 17  -7.138646119543568e-13\n";
    integralFileOneBodyFAD << "      1 38      1 18   2.569912603035684e-13\n";
    integralFileOneBodyFAD << "      1 38      1 19  -2.450935167709958e-13\n";
    integralFileOneBodyFAD << "      1 38      1 20   5.472962024983402e-13\n";
    integralFileOneBodyFAD << "      1 38      1 21   4.259392184660996e-13\n";
    integralFileOneBodyFAD << "      1 38      1 22   7.043464171282988e-13\n";
    integralFileOneBodyFAD << "      1 38      1 23   1.142183379126971e-13\n";
    integralFileOneBodyFAD << "      1 38      1 24  -7.923897192693361e-13\n";
    integralFileOneBodyFAD << "      1 38      1 25   5.663325921504564e-13\n";
    integralFileOneBodyFAD << "      1 38      1 26  -3.640709520967220e-13\n";
    integralFileOneBodyFAD << "      1 38      1 27  -1.832252504016183e-13\n";
    integralFileOneBodyFAD << "      1 38      1 28   1.903638965211618e-13\n";
    integralFileOneBodyFAD << "      1 38      1 29   2.141593835863070e-13\n";
    integralFileOneBodyFAD << "      1 38      1 30  -3.093413318468879e-13\n";
    integralFileOneBodyFAD << "      1 38      1 31  -6.377190533458921e-13\n";
    integralFileOneBodyFAD << "      1 38      1 32  -4.402165107051867e-13\n";
    integralFileOneBodyFAD << "      1 38      1 33  -3.854868904553527e-13\n";
    integralFileOneBodyFAD << "      1 38      1 34   8.685352778778008e-13\n";
    integralFileOneBodyFAD << "      1 38      1 35   1.182635707137718e-13\n";
    integralFileOneBodyFAD << "      1 38      1 36  -3.593118546836929e-13\n";
    integralFileOneBodyFAD << "      1 38      1 37   2.676992294828838e-13\n";
    integralFileOneBodyFAD << "      1 38      1 38   3.133897924626650e+03\n";
    integralFileOneBodyFAD.close();
    //
    integralsOneBodyFAD = MaquisIntegralType{
        {{1, 0, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1}, -2.359242429009664e+03},
        {{1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1}, -2.358797386438359e+03},
        {{1, 2, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1}, -1.444541233850135e+03},
        {{1, 3, 1, 3, -1, -1, -1, -1, -1, -1, -1, -1}, -1.437160009122911e+03},
        {{1, 4, 1, 4, -1, -1, -1, -1, -1, -1, -1, -1}, -6.657064784181758e+03},
        {{1, 5, 1, 5, -1, -1, -1, -1, -1, -1, -1, -1}, -6.128584661601011e+03},
        {{1, 6, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1}, -3.441604675890336e+03},
        {{1, 7, 1, 7, -1, -1, -1, -1, -1, -1, -1, -1}, 1.538026675462226e+03},
        {{1, 8, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1}, 6.026905892315631e+03},
        {{1, 9, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1}, 9.519079627438325e+03},
        {{1, 10, 1, 10, -1, -1, -1, -1, -1, -1, -1, -1}, 1.408367346423375e+03},
        {{1, 11, 1, 11, -1, -1, -1, -1, -1, -1, -1, -1}, 1.877962833530346e+03},
        {{1, 12, 1, 12, -1, -1, -1, -1, -1, -1, -1, -1}, 2.406532356803376e+03},
        {{1, 13, 1, 13, -1, -1, -1, -1, -1, -1, -1, -1}, 2.974933802551137e+03},
        {{1, 14, 1, 14, -1, -1, -1, -1, -1, -1, -1, -1}, 3.590607405365762e+03},
        {{1, 15, 1, 15, -1, -1, -1, -1, -1, -1, -1, -1}, 4.250019051366812e+03},
        {{1, 16, 1, 16, -1, -1, -1, -1, -1, -1, -1, -1}, 4.954168295956210e+03},
        {{1, 17, 1, 17, -1, -1, -1, -1, -1, -1, -1, -1}, 5.702284980641777e+03},
        {{1, 18, 1, 18, -1, -1, -1, -1, -1, -1, -1, -1}, 6.494407879568767e+03},
        {{1, 19, 1, 19, -1, -1, -1, -1, -1, -1, -1, -1}, 7.330063654159827e+03},
        {{1, 20, 1, 20, -1, -1, -1, -1, -1, -1, -1, -1}, 8.209556409406479e+03},
        {{1, 21, 1, 21, -1, -1, -1, -1, -1, -1, -1, -1}, 9.132461752263011e+03},
        {{1, 22, 1, 22, -1, -1, -1, -1, -1, -1, -1, -1}, 1.009901501528866e+03},
        {{1, 23, 1, 23, -1, -1, -1, -1, -1, -1, -1, -1}, 1.110869147300428e+03},
        {{1, 24, 1, 24, -1, -1, -1, -1, -1, -1, -1, -1}, 1.216206565203700e+03},
        {{1, 25, 1, 25, -1, -1, -1, -1, -1, -1, -1, -1}, 1.325854103259667e+03},
        {{1, 26, 1, 26, -1, -1, -1, -1, -1, -1, -1, -1}, 1.439864696031767e+03},
        {{1, 27, 1, 27, -1, -1, -1, -1, -1, -1, -1, -1}, 1.558157142548529e+03},
        {{1, 28, 1, 28, -1, -1, -1, -1, -1, -1, -1, -1}, 1.680827752856561e+03},
        {{1, 29, 1, 29, -1, -1, -1, -1, -1, -1, -1, -1}, 1.807804326901167e+03},
        {{1, 30, 1, 30, -1, -1, -1, -1, -1, -1, -1, -1}, 1.939221830929778e+03},
        {{1, 31, 1, 31, -1, -1, -1, -1, -1, -1, -1, -1}, 2.074605464868975e+03},
        {{1, 32, 1, 32, -1, -1, -1, -1, -1, -1, -1, -1}, 2.214735474735023e+03},
        {{1, 33, 1, 33, -1, -1, -1, -1, -1, -1, -1, -1}, 2.359908634757324e+03},
        {{1, 34, 1, 34, -1, -1, -1, -1, -1, -1, -1, -1}, 2.503011884645919e+03},
        {{1, 35, 1, 35, -1, -1, -1, -1, -1, -1, -1, -1}, 2.651668901080508e+03},
        {{1, 36, 1, 36, -1, -1, -1, -1, -1, -1, -1, -1}, 2.861283145084534e+03},
        {{1, 37, 1, 37, -1, -1, -1, -1, -1, -1, -1, -1}, 2.946541753544523e+03},
        {{1, 38, 1, 38, -1, -1, -1, -1, -1, -1, -1, -1},
         3.133897924626650e+03}};
    // For the two-body FCIDUMP, we kept only the constants > 1 cm-1 for
    // convenience
    integralFileTwoBodyFAD.open("integral_file_test_TwoBodyFAD");
    integralFileTwoBodyFAD << "       1 0       1 0  -4.999867107589869e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1  -3.210749553546130e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2   6.996161115711967e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3   1.571443773098094e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4   2.663750744044597e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5   3.884176384536198e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6   5.751299631128344e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7   6.815690651641560e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8   7.537869346760475e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9   1.630472573510268e+04 \n";
    integralFileTwoBodyFAD << "      1 10      1 10   1.637136661978012e+04 \n";
    integralFileTwoBodyFAD << "       2 0       2 0   2.569781187358162e+02 \n";
    integralFileTwoBodyFAD << "       2 1       2 1   7.706035261773208e+02 \n";
    integralFileTwoBodyFAD << "       2 2       2 2   1.282151911867262e+03 \n";
    integralFileTwoBodyFAD << "       2 3       2 3   1.801678060826892e+03 \n";
    integralFileTwoBodyFAD << "       2 4       2 4   2.275909398565573e+03 \n";
    integralFileTwoBodyFAD << "       2 5       2 5   2.856023548446748e+03 \n";
    integralFileTwoBodyFAD << "       2 6       2 6   3.253239986571359e+03 \n";
    integralFileTwoBodyFAD << "       2 7       2 7   3.958563348401186e+03 \n";
    integralFileTwoBodyFAD << "       2 8       2 8   4.735068233350932e+03 \n";
    integralFileTwoBodyFAD << "       2 9       2 9   5.325000494117511e+03 \n";
    integralFileTwoBodyFAD << "      2 10      2 10   6.857702801147996e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 0       2 0   "
                              "-3.284727244213020e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 0       2 1    "
                              "4.080311386590735e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 0       2 2   "
                              "-2.988080500660043e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 0       2 3   "
                              "-1.868705813829901e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 0       2 0    "
                              "2.551930730716510e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 0       2 1   "
                              "-3.303382387233767e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 0       2 2    "
                              "2.433832430226953e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 0       2 3    "
                              "1.593864453718183e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 0       2 1   "
                              "-1.392468260880556e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 0       2 2    "
                              "1.046968311155213e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 0       2 1   "
                              "-2.726237669356156e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 0       2 2    "
                              "2.141076998089610e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 0       2 1    "
                              "9.898910128722843e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 1       2 0    "
                              "4.080311386590735e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 1       2 2   "
                              "-5.768385840711547e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 1       2 3   "
                              "-5.089352375491885e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 1       2 4   "
                              "-7.153825838075533e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 1       2 5    "
                              "1.029584850004034e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 1       2 6   "
                              "-1.784373335694355e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 1       2 0   "
                              "-3.303382387233766e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 1       2 2    "
                              "4.671150201265141e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 1       2 3    "
                              "4.145709149066609e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 1       2 4    "
                              "5.961941199936081e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 1       2 6    "
                              "1.432939737914276e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 1       2 0   "
                              "-1.392468260880556e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 1       2 2    "
                              "1.970386576852508e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 1       2 3    "
                              "1.783532735367434e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 1       2 4    "
                              "2.728398926900657e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 1       2 0   "
                              "-2.726237669356156e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 1       2 2    "
                              "3.860781707625213e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 1       2 3    "
                              "3.644759010542837e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 1       2 0    "
                              "9.898910128722845e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 1       2 2   "
                              "-1.400968750064177e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 1       2 3   "
                              "-1.326051652580063e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 1       2 2   "
                              "-1.292872672988608e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 0   "
                              "-2.988080500660044e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 1   "
                              "-5.768385840711547e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 2    "
                              "2.801142134498308e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 3    "
                              "6.962081759279173e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 4    "
                              "8.193963583378255e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 5    "
                              "3.378344684665776e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 7    "
                              "2.995437604608070e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 8    "
                              "2.305562408928948e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2       2 9    "
                              "1.497558228422014e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 2      2 10   "
                              "-1.520301519373452e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 0    "
                              "2.433832430226953e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 1    "
                              "4.671150201265141e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 2   "
                              "-2.782059829046955e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 3   "
                              "-5.639209206018207e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 4   "
                              "-6.673225609438829e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 5   "
                              "-2.535838957258751e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 7   "
                              "-2.419368543566544e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 8   "
                              "-1.860265805809907e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2       2 9   "
                              "-1.209869728622875e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 2      2 10    "
                              "1.230416139854112e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 2       2 0    "
                              "1.046968311155213e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 2       2 1    "
                              "1.970386576852508e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 2       2 2   "
                              "-1.905519886316613e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 2       2 3   "
                              "-2.380460784189811e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 2       2 4   "
                              "-2.868054581076223e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 2       2 7   "
                              "-1.013408184746776e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 2       2 0    "
                              "2.141076998089611e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 2       2 1    "
                              "3.860781707625213e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 2       2 3   "
                              "-4.668008245825255e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 2       2 4   "
                              "-5.840263653569885e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 2       2 1   "
                              "-1.400968750064177e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 2       2 3    "
                              "1.692724307676117e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 2       2 4    "
                              "2.122050539075629e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 2       2 1   "
                              "-1.292872672988608e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 2       2 3    "
                              "1.563931832604052e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 0   "
                              "-1.868705813829890e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 1   "
                              "-5.089352375491891e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 2    "
                              "6.962081759279173e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 3   "
                              "-2.038986830634734e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 4    "
                              "8.422847458086603e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 5   "
                              "-3.332178495248326e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 6   "
                              "-5.928934645494640e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 7    "
                              "1.296847184813778e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 8   "
                              "-8.317561695237368e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3       2 9    "
                              "4.374270385581540e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 3      2 10    "
                              "2.349088739269432e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 0    "
                              "1.593864453718169e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 1    "
                              "4.145709149066609e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 2   "
                              "-5.639209206018207e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 3    "
                              "1.576726780865235e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 4   "
                              "-6.823921669935858e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 5    "
                              "2.728030657099655e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 6    "
                              "4.849101303691406e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 7   "
                              "-1.048132260972923e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 8    "
                              "6.745370134377796e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3       2 9   "
                              "-3.536454451630775e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 3      2 10   "
                              "-1.916040882777747e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 1    "
                              "1.783532735367431e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 2   "
                              "-2.380460784189811e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 3    "
                              "5.595564612473115e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 4   "
                              "-2.882497074847803e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 5    "
                              "1.193925048101820e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 6    "
                              "2.106563414627802e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 7   "
                              "-4.387073529628116e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 8    "
                              "2.862872488740680e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 3       2 9   "
                              "-1.482495648159437e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 3       2 1    "
                              "3.644759010542833e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 3       2 2   "
                              "-4.668008245825255e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 3       2 4   "
                              "-5.658183131020966e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 3       2 5    "
                              "2.541196880784551e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 3       2 6    "
                              "4.303859263961091e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 3       2 1   "
                              "-1.326051652580063e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 3       2 2    "
                              "1.692724307676117e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 3       2 4    "
                              "2.051003749055949e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 3       2 6   "
                              "-1.537535523632375e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 3       2 2    "
                              "1.563931832604052e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 3       2 4    "
                              "1.896344383368165e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 1   "
                              "-7.153825838075504e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 2    "
                              "8.193963583378255e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 3    "
                              "8.422847458086603e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 4    "
                              "1.417012213342013e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 5   "
                              "-7.298487898570368e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 6   "
                              "-2.831449731364563e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 7   "
                              "-3.566075502287342e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 8   "
                              "-3.549594182609876e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4       2 9   "
                              "-1.332877396871099e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 4      2 10    "
                              "1.597413108566536e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 1    "
                              "5.961941199936073e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 2   "
                              "-6.673225609438829e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 3   "
                              "-6.823921669935858e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 4   "
                              "-1.155886488195299e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 5    "
                              "5.913558666520155e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 6    "
                              "2.301742658267990e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 7    "
                              "2.870199137348181e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 8    "
                              "2.865548953353978e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4       2 9    "
                              "1.077288158496926e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 4      2 10   "
                              "-1.295549003315189e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 1    "
                              "2.728398926900647e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 2   "
                              "-2.868054581076223e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 3   "
                              "-2.882497074847803e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 4   "
                              "-4.993945515708544e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 5    "
                              "2.498254853469479e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 6    "
                              "9.827606592886538e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 7    "
                              "1.188386623087850e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 8    "
                              "1.198793496548528e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4       2 9    "
                              "4.517805884059679e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 4      2 10   "
                              "-5.492070086242443e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 4       2 2   "
                              "-5.840263653569885e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 4       2 3   "
                              "-5.658183131020966e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 4       2 4   "
                              "-1.027814612081335e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 4       2 5    "
                              "4.899019466432280e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 4       2 6    "
                              "1.967239206932401e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 4       2 7    "
                              "2.269065133179009e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 4       2 8    "
                              "2.336647615369945e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 4      2 10   "
                              "-1.086434005727392e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 4       2 2    "
                              "2.122050539075629e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 4       2 3    "
                              "2.051003749055949e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 4       2 4    "
                              "3.737094920055676e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 4       2 5   "
                              "-1.773536430507336e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 4       2 6   "
                              "-7.119597352753580e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 4       2 3    "
                              "1.896344383368165e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 4       2 5   "
                              "-1.642475989091999e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 1    "
                              "1.029584850004006e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 2    "
                              "3.378344684665812e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 3   "
                              "-3.332178495248316e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 4   "
                              "-7.298487898570368e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 5   "
                              "-6.703343363960392e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 6    "
                              "8.894887145567910e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 7    "
                              "1.703981792409626e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 8   "
                              "-1.494338497685526e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5       2 9    "
                              "4.017765166589879e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 5      2 10   "
                              "-6.401027202004260e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5       2 2   "
                              "-2.535838957258765e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5       2 3    "
                              "2.728030657099658e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5       2 4    "
                              "5.913558666520154e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5       2 5    "
                              "5.419833292430021e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5       2 6   "
                              "-7.211522100071874e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5       2 7   "
                              "-1.380397524981230e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5       2 8    "
                              "1.250256146411232e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5       2 9   "
                              "-3.250758670644857e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 5      2 10    "
                              "5.151207797076097e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 5       2 3    "
                              "1.193925048101820e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 5       2 4    "
                              "2.498254853469479e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 5       2 5    "
                              "2.273767402683972e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 5       2 6   "
                              "-3.052896196016509e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 5       2 7   "
                              "-5.822809404803636e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 5       2 8    "
                              "5.777886299654141e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 5       2 9   "
                              "-1.366111317915933e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 5      2 10    "
                              "2.130802009503953e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 5       2 3    "
                              "2.541196880784548e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 5       2 4    "
                              "4.899019466432280e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 5       2 5    "
                              "4.395434511131093e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 5       2 6   "
                              "-6.012988730453159e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 5       2 7   "
                              "-1.130060948089803e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 5       2 8    "
                              "1.268994946619731e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 5       2 9   "
                              "-2.641232153392851e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 5       2 4   "
                              "-1.773536430507336e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 5       2 5   "
                              "-1.591035244375780e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 5       2 6    "
                              "2.177318971603242e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 5       2 7    "
                              "4.054464594139078e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 5       2 4   "
                              "-1.642475989091999e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 5       2 5   "
                              "-1.468996738069974e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 5       2 6    "
                              "2.017684486939142e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6       2 1   "
                              "-1.784373335694347e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6       2 3   "
                              "-5.928934645494635e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6       2 4   "
                              "-2.831449731364564e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6       2 5    "
                              "8.894887145567910e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6       2 6    "
                              "1.197467014438249e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6       2 7   "
                              "-3.853378354587493e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6       2 8   "
                              "-3.571735584404719e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6       2 9   "
                              "-4.494136971784449e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 6      2 10    "
                              "6.499554000425073e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6       2 1    "
                              "1.432939737914270e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6       2 3    "
                              "4.849101303691411e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6       2 4    "
                              "2.301742658267989e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6       2 5   "
                              "-7.211522100071874e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6       2 6   "
                              "-9.718822990700186e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6       2 7    "
                              "3.123201294521006e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6       2 8    "
                              "2.908005904154444e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6       2 9    "
                              "3.631155807796803e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 6      2 10   "
                              "-5.280869161006557e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 6       2 3    "
                              "2.106563414627800e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 6       2 4    "
                              "9.827606592886535e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 6       2 5   "
                              "-3.052896196016509e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 6       2 6   "
                              "-4.128850191326192e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 6       2 7    "
                              "1.320567810911419e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 6       2 8    "
                              "1.246636358465198e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 6       2 9    "
                              "1.520722898415504e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 6      2 10   "
                              "-2.250794069603901e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 6       2 3    "
                              "4.303859263961091e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 6       2 4    "
                              "1.967239206932402e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 6       2 5   "
                              "-6.012988730453159e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 6       2 6   "
                              "-8.192923728174321e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 6       2 7    "
                              "2.590026409798208e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 6       2 8    "
                              "2.501286650343352e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 6       2 9    "
                              "2.941297609416032e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 6      2 10   "
                              "-4.487282447436260e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 6       2 3   "
                              "-1.537535523632375e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 6       2 4   "
                              "-7.119597352753580e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 6       2 5    "
                              "2.177318971603242e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 6       2 6    "
                              "2.968128860639846e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 6       2 7   "
                              "-9.359069647205242e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 6       2 8   "
                              "-9.005928584056226e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 6       2 9   "
                              "-1.070248767816968e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 6      2 10    "
                              "1.622634053486745e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 6       2 5    "
                              "2.017684486939142e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 6       2 6    "
                              "2.753394733772488e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7       2 2    "
                              "2.995437604608098e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7       2 3    "
                              "1.296847184813781e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7       2 4   "
                              "-3.566075502287342e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7       2 5    "
                              "1.703981792409626e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7       2 6   "
                              "-3.853378354587493e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7       2 7   "
                              "-1.526496535878683e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7       2 8    "
                              "2.372456045292365e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7       2 9    "
                              "1.756129097567093e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 7      2 10   "
                              "-3.200356074943027e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7       2 2   "
                              "-2.419368543566544e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7       2 3   "
                              "-1.048132260972921e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7       2 4    "
                              "2.870199137348181e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7       2 5   "
                              "-1.380397524981229e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7       2 6    "
                              "3.123201294521006e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7       2 7    "
                              "1.235901989280693e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7       2 8   "
                              "-1.925829476233712e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7       2 9   "
                              "-1.426049970748351e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 7      2 10    "
                              "2.592367564713680e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7       2 2   "
                              "-1.013408184746767e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7       2 3   "
                              "-4.387073529628124e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7       2 4    "
                              "1.188386623087850e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7       2 5   "
                              "-5.822809404803636e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7       2 6    "
                              "1.320567810911419e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7       2 7    "
                              "5.207536126860404e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7       2 8   "
                              "-8.182168002205952e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7       2 9   "
                              "-6.061114169503517e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 7      2 10    "
                              "1.094561935310272e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 7       2 4    "
                              "2.269065133179009e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 7       2 5   "
                              "-1.130060948089804e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 7       2 6    "
                              "2.590026409798208e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 7       2 7    "
                              "1.014494595332452e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 7       2 8   "
                              "-1.619131998550263e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 7       2 9   "
                              "-1.191198100873517e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 7      2 10    "
                              "2.147959907634376e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 7       2 5    "
                              "4.054464594139078e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 7       2 6   "
                              "-9.359069647205242e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 7       2 7   "
                              "-3.666759324278461e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 7       2 8    "
                              "5.848176844070152e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 7       2 9    "
                              "4.253366829643912e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 7       2 7   "
                              "-3.396837755305315e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8       2 2    "
                              "2.305562408928947e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8       2 3   "
                              "-8.317561695237353e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8       2 4   "
                              "-3.549594182609876e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8       2 5   "
                              "-1.494338497685525e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8       2 6   "
                              "-3.571735584404719e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8       2 7    "
                              "2.372456045292366e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8       2 8    "
                              "2.355734664584248e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8       2 9   "
                              "-1.858796588023494e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 8      2 10    "
                              "2.735145396463816e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8       2 2   "
                              "-1.860265805809907e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8       2 3    "
                              "6.745370134377781e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8       2 4    "
                              "2.865548953353978e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8       2 5    "
                              "1.250256146411235e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8       2 6    "
                              "2.908005904154444e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8       2 7   "
                              "-1.925829476233712e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8       2 8   "
                              "-1.913770103311484e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8       2 9    "
                              "1.508134375274986e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 8      2 10   "
                              "-2.234716112660251e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 8       2 3    "
                              "2.862872488740680e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 8       2 4    "
                              "1.198793496548527e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 8       2 5    "
                              "5.777886299654140e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 8       2 6    "
                              "1.246636358465198e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 8       2 7   "
                              "-8.182168002205952e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 8       2 8   "
                              "-8.152060972764459e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 8       2 9    "
                              "6.397272639737099e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 8      2 10   "
                              "-9.673090992867571e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 8       2 4    "
                              "2.336647615369939e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 8       2 5    "
                              "1.268994946619731e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 8       2 6    "
                              "2.501286650343352e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 8       2 7   "
                              "-1.619131998550264e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 8       2 8   "
                              "-1.622133415357714e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 8       2 9    "
                              "1.261371567272619e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 8      2 10   "
                              "-1.959826459762753e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 8       2 6   "
                              "-9.005928584056230e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 8       2 7    "
                              "5.848176844070152e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 8       2 8    "
                              "5.862926376735706e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 8       2 9   "
                              "-4.552211551146165e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 8      2 10    "
                              "7.001489898205741e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 8       2 8    "
                              "5.463276728836961e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9       2 2    "
                              "1.497558228422007e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9       2 3    "
                              "4.374270385581537e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9       2 4   "
                              "-1.332877396871099e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9       2 5    "
                              "4.017765166589879e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9       2 6   "
                              "-4.494136971784449e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9       2 7    "
                              "1.756129097567093e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9       2 8   "
                              "-1.858796588023494e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9       2 9   "
                              "-1.924757249600509e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 0       2 9      2 10   "
                              "-7.965785749009591e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9       2 2   "
                              "-1.209869728622874e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9       2 3   "
                              "-3.536454451630782e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9       2 4    "
                              "1.077288158496927e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9       2 5   "
                              "-3.250758670644853e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9       2 6    "
                              "3.631155807796803e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9       2 7   "
                              "-1.426049970748351e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9       2 8    "
                              "1.508134375274986e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9       2 9    "
                              "1.559282784835807e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 2       2 9      2 10    "
                              "6.479802931818537e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 9       2 3   "
                              "-1.482495648159437e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 9       2 4    "
                              "4.517805884059674e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 9       2 5   "
                              "-1.366111317915933e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 9       2 6    "
                              "1.520722898415503e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 9       2 7   "
                              "-6.061114169503514e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 9       2 8    "
                              "6.397272639737099e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 9       2 9    "
                              "6.581877662889448e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 4       2 9      2 10    "
                              "2.769873474544293e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 9       2 5   "
                              "-2.641232153392851e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 9       2 6    "
                              "2.941297609416032e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 9       2 7   "
                              "-1.191198100873517e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 9       2 8    "
                              "1.261371567272619e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 9       2 9    "
                              "1.284858463010085e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 6       2 9      2 10    "
                              "5.517456553623747e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 9       2 6   "
                              "-1.070248767816967e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 9       2 7    "
                              "4.253366829643913e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 9       2 8   "
                              "-4.552211551146165e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 9       2 9   "
                              "-4.633511357968788e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 8       2 9      2 10   "
                              "-1.979024442684471e+00 \n";
    integralFileTwoBodyFAD << "       1 0      1 10       2 9       2 9   "
                              "-4.303707474426987e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10       2 2   "
                              "-1.520301519373445e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10       2 3    "
                              "2.349088739269432e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10       2 4    "
                              "1.597413108566538e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10       2 5   "
                              "-6.401027202004261e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10       2 6    "
                              "6.499554000425071e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10       2 7   "
                              "-3.200356074943029e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10       2 8    "
                              "2.735145396463815e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10       2 9   "
                              "-7.965785749009592e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 0      2 10      2 10    "
                              "3.103690442760890e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10       2 2    "
                              "1.230416139854110e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10       2 3   "
                              "-1.916040882777747e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10       2 4   "
                              "-1.295549003315189e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10       2 5    "
                              "5.151207797076102e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10       2 6   "
                              "-5.280869161006556e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10       2 7    "
                              "2.592367564713681e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10       2 8   "
                              "-2.234716112660252e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10       2 9    "
                              "6.479802931818537e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 2      2 10      2 10   "
                              "-2.525337354818048e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 4      2 10       2 4   "
                              "-5.492070086242442e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4      2 10       2 5    "
                              "2.130802009503949e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 4      2 10       2 6   "
                              "-2.250794069603901e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4      2 10       2 7    "
                              "1.094561935310272e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4      2 10       2 8   "
                              "-9.673090992867571e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4      2 10       2 9    "
                              "2.769873474544295e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 4      2 10      2 10   "
                              "-1.080485149673152e+03 \n";
    integralFileTwoBodyFAD << "       1 0       1 6      2 10       2 4   "
                              "-1.086434005727393e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6      2 10       2 6   "
                              "-4.487282447436260e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6      2 10       2 7    "
                              "2.147959907634375e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6      2 10       2 8   "
                              "-1.959826459762754e+01 \n";
    integralFileTwoBodyFAD << "       1 0       1 6      2 10       2 9    "
                              "5.517456553623748e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 6      2 10      2 10   "
                              "-2.160981958127234e+02 \n";
    integralFileTwoBodyFAD << "       1 0       1 8      2 10       2 6    "
                              "1.622634053486745e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8      2 10       2 8    "
                              "7.001489898205741e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8      2 10       2 9   "
                              "-1.979024442684471e+00 \n";
    integralFileTwoBodyFAD << "       1 0       1 8      2 10      2 10    "
                              "7.785748816323367e+01 \n";
    integralFileTwoBodyFAD << "       1 0      1 10      2 10      2 10    "
                              "7.304753787658052e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 0       2 0   "
                              "-4.457629551969945e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 0       2 1    "
                              "5.586879549865388e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 0       2 2   "
                              "-4.096836643331162e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 0       2 3   "
                              "-2.589416118903074e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 0       2 0   "
                              "-2.633211657127593e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 0       2 1    "
                              "3.544633365081027e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 0       2 2   "
                              "-2.625956126994469e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 0       2 3   "
                              "-1.781442696689945e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 0       2 1    "
                              "8.823396009909217e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 0       2 2   "
                              "-6.743992800579903e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 0       2 1   "
                              "-2.895320936801980e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 0       2 2    "
                              "2.274405728380943e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 0       2 1   "
                              "-1.730807465470138e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 1       2 0    "
                              "5.586879549865388e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 1       2 2   "
                              "-7.898666115135379e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 1       2 3   "
                              "-6.977957227511713e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 1       2 4   "
                              "-9.859833622090935e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 1       2 5    "
                              "1.415339823194447e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 1       2 6   "
                              "-2.438705175671903e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 1       2 8   "
                              "-1.072460786410635e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 1       2 0    "
                              "3.544633365081026e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 1       2 2   "
                              "-5.013246559842777e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 1       2 3   "
                              "-4.473105088136514e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 1       2 4   "
                              "-6.547100969860350e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 1       2 6   "
                              "-1.528354157962007e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 1       2 0    "
                              "8.823396009909216e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 1       2 2   "
                              "-1.249099522197051e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 1       2 3   "
                              "-1.148745155732089e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 1       2 4   "
                              "-1.822899804609057e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 1       2 0   "
                              "-2.895320936801980e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 1       2 2    "
                              "4.099039183720793e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 1       2 3    "
                              "3.870940551514999e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 1       2 0   "
                              "-1.730807465470138e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 1       2 2    "
                              "2.451597067677716e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 0   "
                              "-4.096836643331162e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 1   "
                              "-7.898666115135379e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 2    "
                              "4.026809655968736e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 3    "
                              "9.533740049829545e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 4    "
                              "1.123415262635888e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 5    "
                              "4.550102049286187e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 7    "
                              "4.099236223256630e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 8    "
                              "3.154459895987023e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2       2 9    "
                              "2.049532033323135e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 2      2 10   "
                              "-2.081479198391211e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 0   "
                              "-2.625956126994469e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 1   "
                              "-5.013246559842776e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 2    "
                              "3.486206322710915e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 3    "
                              "6.053411383368234e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 4    "
                              "7.198349927955807e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 5    "
                              "2.546156025363794e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 7    "
                              "2.591506863007371e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 8    "
                              "1.990500127001256e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2       2 9    "
                              "1.296047570331972e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 2      2 10   "
                              "-1.320012570275538e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 2       2 0   "
                              "-6.743992800579902e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 2       2 1   "
                              "-1.249099522197051e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 2       2 2    "
                              "1.588634273425861e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 2       2 3    "
                              "1.509761598759491e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 2       2 4    "
                              "1.845337029158174e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 2       2 0    "
                              "2.274405728380945e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 2       2 1    "
                              "4.099039183720794e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 2       2 3   "
                              "-4.954502197377938e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 2       2 4   "
                              "-6.199420936168715e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 2       2 1    "
                              "2.451597067677716e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 2       2 3   "
                              "-2.964881898903597e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 0   "
                              "-2.589416118903089e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 1   "
                              "-6.977957227511712e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 2    "
                              "9.533740049829545e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 3   "
                              "-2.764271921521222e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 4    "
                              "1.153465525858069e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 5   "
                              "-4.573771834539674e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 6   "
                              "-8.136675743511871e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 7    "
                              "1.775043874035217e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 8   "
                              "-1.139279902052941e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3       2 9    "
                              "5.987603655397630e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 3      2 10    "
                              "3.221803047775054e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 0   "
                              "-1.781442696689945e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 1   "
                              "-4.473105088136510e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 2    "
                              "6.053411383368235e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 3   "
                              "-1.620015215629235e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 4    "
                              "7.326486073470033e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 5   "
                              "-2.957236128215969e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 6   "
                              "-5.246678329675051e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 7    "
                              "1.122598623169607e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 8   "
                              "-7.251075798011815e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3       2 9    "
                              "3.789217507298867e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 3      2 10    "
                              "2.069022110495699e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 3       2 1   "
                              "-1.148745155732090e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 3       2 2    "
                              "1.509761598759491e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 3       2 3   "
                              "-3.000836027715816e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 3       2 4    "
                              "1.829035662226728e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 3       2 5   "
                              "-7.801003823052897e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 3       2 6   "
                              "-1.362561801522246e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 3       2 7    "
                              "2.760620515744176e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 3       2 8   "
                              "-1.825960811296100e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 3       2 1    "
                              "3.870940551514999e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 3       2 2   "
                              "-4.954502197377938e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 3       2 4   "
                              "-6.004325429104954e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 3       2 5    "
                              "2.708858206606646e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 3       2 6    "
                              "4.534297069014713e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 3       2 2   "
                              "-2.964881898903597e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 3       2 4   "
                              "-3.593929799520220e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 1   "
                              "-9.859833622090934e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 2    "
                              "1.123415262635889e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 3    "
                              "1.153465525858069e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 4    "
                              "1.943457186380929e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 5   "
                              "-9.995132158315340e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 6   "
                              "-3.880453144214783e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 7   "
                              "-4.876384541425468e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 8   "
                              "-4.857049433989589e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4       2 9   "
                              "-1.824334704467028e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 4      2 10    "
                              "2.188087307654219e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 1   "
                              "-6.547100969860328e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 2    "
                              "7.198349927955810e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 3    "
                              "7.326486073470033e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 4    "
                              "1.248631492851935e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 5   "
                              "-6.349332585316281e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 6   "
                              "-2.478479621974622e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 7   "
                              "-3.064871319823028e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 8   "
                              "-3.068365192265271e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4       2 9   "
                              "-1.154350779134183e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 4      2 10    "
                              "1.392298577704565e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 1   "
                              "-1.822899804609048e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 2    "
                              "1.845337029158173e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 3    "
                              "1.829035662226727e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 4    "
                              "3.226385831330300e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 5   "
                              "-1.585052114112287e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 6   "
                              "-6.286568347067821e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 7   "
                              "-7.436477558802723e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 8   "
                              "-7.564426310747179e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4       2 9   "
                              "-2.852205906660336e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 4      2 10    "
                              "3.495015082827121e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 4       2 2   "
                              "-6.199420936168717e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 4       2 3   "
                              "-6.004325429104954e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 4       2 4   "
                              "-1.091113269607886e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 4       2 5    "
                              "5.195899699454709e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 4       2 6    "
                              "2.085247865990196e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 4       2 7    "
                              "2.426654210686763e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 4       2 8    "
                              "2.492485166653142e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 4      2 10   "
                              "-1.152958553230018e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 4       2 3   "
                              "-3.593929799520220e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 4       2 5    "
                              "3.114070807786512e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 4       2 6    "
                              "1.247838238673187e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 1    "
                              "1.415339823194463e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 2    "
                              "4.550102049286243e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 3   "
                              "-4.573771834539691e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 4   "
                              "-9.995132158315341e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 5   "
                              "-9.175786048075314e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 6    "
                              "1.218304432612649e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 7    "
                              "2.333501145760420e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 8   "
                              "-2.061461736284479e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5       2 9    "
                              "5.500561753150241e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 5      2 10   "
                              "-8.752766262950027e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5       2 2    "
                              "2.546156025363766e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5       2 3   "
                              "-2.957236128215967e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5       2 4   "
                              "-6.349332585316281e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5       2 5   "
                              "-5.808308224762369e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5       2 6    "
                              "7.747258227388838e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5       2 7    "
                              "1.481572049056897e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5       2 8   "
                              "-1.376902317603778e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5       2 9    "
                              "3.485434155131459e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 5      2 10   "
                              "-5.499337648114057e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 5       2 3   "
                              "-7.801003823052904e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 5       2 4   "
                              "-1.585052114112287e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 5       2 5   "
                              "-1.434643881367362e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 5       2 6    "
                              "1.940177984279219e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 5       2 7    "
                              "3.685733728801262e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 5       2 8   "
                              "-3.881736394290436e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 5       2 9    "
                              "8.626285889977360e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 5      2 10   "
                              "-1.331368796000010e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 5       2 3    "
                              "2.708858206606646e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 5       2 4    "
                              "5.195899699454709e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 5       2 5    "
                              "4.662981931487190e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 5       2 6   "
                              "-6.377433945415035e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 5       2 7   "
                              "-1.194218870577143e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 5       2 8    "
                              "1.316093853223485e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 5       2 9   "
                              "-2.799012360995751e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 5       2 4    "
                              "3.114070807786512e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 5       2 5    "
                              "2.799026088724461e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 5       2 6   "
                              "-3.819670485141559e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6       2 1   "
                              "-2.438705175671902e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6       2 3   "
                              "-8.136675743511879e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6       2 4   "
                              "-3.880453144214783e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6       2 5    "
                              "1.218304432612649e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6       2 6    "
                              "1.640519612651926e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6       2 7   "
                              "-5.277519146452846e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6       2 8   "
                              "-4.896726243101808e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6       2 9   "
                              "-6.150782574308525e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 6      2 10    "
                              "8.906504831355900e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6       2 1   "
                              "-1.528354157962000e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6       2 3   "
                              "-5.246678329675051e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6       2 4   "
                              "-2.478479621974622e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6       2 5    "
                              "7.747258227388839e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6       2 6    "
                              "1.045076850443367e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6       2 7   "
                              "-3.354159511611819e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6       2 8   "
                              "-3.134849417383397e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6       2 9   "
                              "-3.889522473432190e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 6      2 10    "
                              "5.683666399448378e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 6       2 3   "
                              "-1.362561801522247e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 6       2 4   "
                              "-6.286568347067818e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 6       2 5    "
                              "1.940177984279219e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 6       2 6    "
                              "2.631392050094647e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 6       2 7   "
                              "-8.382115338489656e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 6       2 8   "
                              "-7.991987013406877e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 6       2 9   "
                              "-9.587835720001385e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 6      2 10    "
                              "1.437706639836913e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 6       2 3    "
                              "4.534297069014717e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 6       2 4    "
                              "2.085247865990196e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 6       2 5   "
                              "-6.377433945415035e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 6       2 6   "
                              "-8.689957917015739e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 6       2 7    "
                              "2.744790145473134e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 6       2 8    "
                              "2.645144412229412e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 6       2 9    "
                              "3.127256675103653e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 6      2 10   "
                              "-4.755186837855253e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 6       2 4    "
                              "1.247838238673187e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 6       2 5   "
                              "-3.819670485141559e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 6       2 6   "
                              "-5.198931969245064e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 6       2 7    "
                              "1.647355097856701e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 6       2 8    "
                              "1.591502419242160e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7       2 2    "
                              "4.099236223256686e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7       2 3    "
                              "1.775043874035209e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7       2 4   "
                              "-4.876384541425474e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7       2 5    "
                              "2.333501145760420e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7       2 6   "
                              "-5.277519146452846e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7       2 7   "
                              "-2.090160458286823e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7       2 8    "
                              "3.250376078222909e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7       2 9    "
                              "2.406193716203966e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 7      2 10   "
                              "-4.382543647821165e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7       2 2    "
                              "2.591506863007371e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7       2 3    "
                              "1.122598623169606e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7       2 4   "
                              "-3.064871319823024e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7       2 5    "
                              "1.481572049056897e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7       2 6   "
                              "-3.354159511611819e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7       2 7   "
                              "-1.326046253590376e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7       2 8    "
                              "2.070948914317344e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7       2 9    "
                              "1.533730031004696e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 7      2 10   "
                              "-2.782946462827163e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 7       2 3    "
                              "2.760620515744186e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 7       2 4   "
                              "-7.436477558802723e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 7       2 5    "
                              "3.685733728801262e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 7       2 6   "
                              "-8.382115338489656e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 7       2 7   "
                              "-3.296506992663207e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 7       2 8    "
                              "5.212488665443411e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 7       2 9    "
                              "3.858025883183865e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 7      2 10   "
                              "-6.943714719025521e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 7       2 4    "
                              "2.426654210686765e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 7       2 5   "
                              "-1.194218870577142e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 7       2 6    "
                              "2.744790145473134e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 7       2 7    "
                              "1.075387200487358e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 7       2 8   "
                              "-1.715197084999770e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 7       2 9   "
                              "-1.256027146080472e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 7      2 10    "
                              "2.279602297199015e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 7       2 6    "
                              "1.647355097856701e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 7       2 7    "
                              "6.456504337567208e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 7       2 8   "
                              "-1.029351659267909e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 1   "
                              "-1.072460786410652e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 2    "
                              "3.154459895987020e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 3   "
                              "-1.139279902052948e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 4   "
                              "-4.857049433989582e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 5   "
                              "-2.061461736284479e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 6   "
                              "-4.896726243101812e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 7    "
                              "3.250376078222909e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 8    "
                              "3.228032631788115e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8       2 9   "
                              "-2.546363979980296e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 8      2 10    "
                              "3.752790329065256e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8       2 2    "
                              "1.990500127001258e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8       2 3   "
                              "-7.251075798011831e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8       2 4   "
                              "-3.068365192265271e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8       2 5   "
                              "-1.376902317603778e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8       2 6   "
                              "-3.134849417383397e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8       2 7    "
                              "2.070948914317344e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8       2 8    "
                              "2.059427107508353e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8       2 9   "
                              "-1.621077623892800e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 8      2 10    "
                              "2.415562839422608e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 8       2 3   "
                              "-1.825960811296103e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 8       2 4   "
                              "-7.564426310747179e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 8       2 5   "
                              "-3.881736394290450e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 8       2 6   "
                              "-7.991987013406877e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 8       2 7    "
                              "5.212488665443411e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 8       2 8    "
                              "5.204120180804814e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 8       2 9   "
                              "-4.070085339804571e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 8      2 10    "
                              "6.237789957088010e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 8       2 4    "
                              "2.492485166653144e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 8       2 5    "
                              "1.316093853223485e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 8       2 6    "
                              "2.645144412229412e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 8       2 7   "
                              "-1.715197084999770e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 8       2 8   "
                              "-1.718651925716283e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 8       2 9    "
                              "1.335862154989681e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 8      2 10   "
                              "-2.064903844978512e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 8       2 6    "
                              "1.591502419242161e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 8       2 7   "
                              "-1.029351659267909e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 8       2 8   "
                              "-1.030555941314444e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 8      2 10   "
                              "-1.252907494784046e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9       2 2    "
                              "2.049532033323142e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9       2 3    "
                              "5.987603655397601e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9       2 4   "
                              "-1.824334704467028e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9       2 5    "
                              "5.500561753150241e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9       2 6   "
                              "-6.150782574308521e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9       2 7    "
                              "2.406193716203965e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9       2 8   "
                              "-2.546363979980296e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9       2 9   "
                              "-2.635835273857021e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 1       2 9      2 10   "
                              "-1.091867807323455e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9       2 2    "
                              "1.296047570331955e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9       2 3    "
                              "3.789217507298867e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9       2 4   "
                              "-1.154350779134183e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9       2 5    "
                              "3.485434155131459e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9       2 6   "
                              "-3.889522473432191e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9       2 7    "
                              "1.533730031004695e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9       2 8   "
                              "-1.621077623892800e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9       2 9   "
                              "-1.673834988055630e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 3       2 9      2 10   "
                              "-6.979851870801234e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 9       2 4   "
                              "-2.852205906660336e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 9       2 5    "
                              "8.626285889977350e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 9       2 6   "
                              "-9.587835720001388e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 9       2 7    "
                              "3.858025883183865e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 9       2 8   "
                              "-4.070085339804571e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 9       2 9   "
                              "-4.171402579764399e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 5       2 9      2 10   "
                              "-1.771431285625149e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 9       2 5   "
                              "-2.799012360995746e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 9       2 6    "
                              "3.127256675103653e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 9       2 7   "
                              "-1.256027146080472e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 9       2 8    "
                              "1.335862154989681e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 9       2 9    "
                              "1.360631752145197e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 7       2 9      2 10    "
                              "5.826741370141407e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9       2 9       2 9    "
                              "8.178121829989706e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10       2 2   "
                              "-2.081479198391208e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10       2 3    "
                              "3.221803047775058e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10       2 4    "
                              "2.188087307654217e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10       2 5   "
                              "-8.752766262950027e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10       2 6    "
                              "8.906504831355915e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10       2 7   "
                              "-4.382543647821165e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10       2 8    "
                              "3.752790329065256e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10       2 9   "
                              "-1.091867807323455e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 1      2 10      2 10    "
                              "4.254449857378064e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10       2 2   "
                              "-1.320012570275538e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10       2 3    "
                              "2.069022110495700e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10       2 4    "
                              "1.392298577704564e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10       2 5   "
                              "-5.499337648114057e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10       2 6    "
                              "5.683666399448375e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10       2 7   "
                              "-2.782946462827164e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10       2 8    "
                              "2.415562839422608e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10       2 9   "
                              "-6.979851870801234e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 3      2 10      2 10    "
                              "2.720876018621629e+03 \n";
    integralFileTwoBodyFAD << "       1 1       1 5      2 10       2 4    "
                              "3.495015082827121e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5      2 10       2 5   "
                              "-1.331368796000010e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5      2 10       2 6    "
                              "1.437706639836913e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5      2 10       2 7   "
                              "-6.943714719025521e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 5      2 10       2 8    "
                              "6.237789957088010e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5      2 10       2 9   "
                              "-1.771431285625149e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 5      2 10      2 10    "
                              "6.917099271935310e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 7      2 10       2 4   "
                              "-1.152958553230019e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7      2 10       2 6   "
                              "-4.755186837855250e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7      2 10       2 7    "
                              "2.279602297199016e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7      2 10       2 8   "
                              "-2.064903844978512e+01 \n";
    integralFileTwoBodyFAD << "       1 1       1 7      2 10       2 9    "
                              "5.826741370141407e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 7      2 10      2 10   "
                              "-2.286064961767198e+02 \n";
    integralFileTwoBodyFAD << "       1 1       1 9      2 10       2 8   "
                              "-1.252907494784047e+00 \n";
    integralFileTwoBodyFAD << "       1 1       1 9      2 10      2 10   "
                              "-1.375296255861435e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 0       2 0    "
                              "2.551930730716524e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 0       2 1   "
                              "-3.303382387233766e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 0       2 2    "
                              "2.433832430226953e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 0       2 3    "
                              "1.593864453718176e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 0       2 0   "
                              "-3.654733389030808e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 0       2 1    "
                              "4.735919841079972e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 0       2 2   "
                              "-3.489783038973312e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 0       2 3   "
                              "-2.277646296639988e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 0       2 0   "
                              "-3.012580490735540e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 0       2 1    "
                              "4.096365102835440e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 0       2 2   "
                              "-3.038853937095198e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 0       2 3   "
                              "-2.065480010246626e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 0       2 1    "
                              "9.847171376405946e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 0       2 2   "
                              "-7.571597286196542e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 0       2 1   "
                              "-3.151543766714967e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 0       2 2    "
                              "2.481706819616804e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 0       2 1   "
                              "-3.125876669004097e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 1       2 0   "
                              "-3.303382387233766e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 1       2 2    "
                              "4.671150201265141e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 1       2 3    "
                              "4.145709149066609e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 1       2 4    "
                              "5.961941199936081e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 1       2 6    "
                              "1.432939737914305e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 1       2 0    "
                              "4.735919841079972e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 1       2 2   "
                              "-6.696693101508555e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 1       2 3   "
                              "-5.944120868937016e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 1       2 4   "
                              "-8.531917225602516e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 1       2 5    "
                              "1.211427258319006e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 1       2 6   "
                              "-2.056720199563080e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 1       2 0    "
                              "4.096365102835439e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 1       2 2   "
                              "-5.793608983847445e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 1       2 3   "
                              "-5.176184784101910e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 1       2 4   "
                              "-7.580613009423137e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 1       2 5    "
                              "1.060793256915073e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 1       2 6   "
                              "-1.767087056018776e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 1       2 0    "
                              "9.847171376405947e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 1       2 2   "
                              "-1.394109019314518e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 1       2 3   "
                              "-1.289598697485622e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 1       2 4   "
                              "-2.054573374740656e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 1       2 0   "
                              "-3.151543766714967e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 1       2 2    "
                              "4.461175351248311e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 1       2 3    "
                              "4.222429427656757e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 1       2 0   "
                              "-3.125876669004097e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 1       2 2    "
                              "4.426763503373093e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 0    "
                              "2.433832430226953e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 1    "
                              "4.671150201265141e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 2   "
                              "-2.782059829046955e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 3   "
                              "-5.639209206018207e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 4   "
                              "-6.673225609438826e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 5   "
                              "-2.535838957258808e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 7   "
                              "-2.419368543566544e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 8   "
                              "-1.860265805809944e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2       2 9   "
                              "-1.209869728622871e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 2      2 10    "
                              "1.230416139854118e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 0   "
                              "-3.489783038973311e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 1   "
                              "-6.696693101508555e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 2    "
                              "4.003871181461451e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 3    "
                              "8.084347622615087e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 4    "
                              "9.567434689579468e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 5    "
                              "3.654065328147420e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 7    "
                              "3.469664676211188e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 8    "
                              "2.667446337745307e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2       2 9    "
                              "1.734837092442033e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 2      2 10   "
                              "-1.764152790431526e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 0   "
                              "-3.038853937095195e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 1   "
                              "-5.793608983847445e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 2    "
                              "4.171563983873237e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 3    "
                              "6.995720722709600e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 4    "
                              "8.328495768830192e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 5    "
                              "2.926922722445232e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 7    "
                              "2.995335256302905e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 8    "
                              "2.299637577010111e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2       2 9    "
                              "1.497673462499823e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 2      2 10   "
                              "-1.525628551227807e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 2       2 0   "
                              "-7.571597286196542e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 2       2 1   "
                              "-1.394109019314518e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 2       2 2    "
                              "1.929094579319154e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 2       2 3    "
                              "1.685115800648999e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 2       2 4    "
                              "2.070443939231003e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 2       2 0    "
                              "2.481706819616804e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 2       2 1    "
                              "4.461175351248311e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 2       2 3   "
                              "-5.391394663410128e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 2       2 4   "
                              "-6.758362232961114e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 2       2 1    "
                              "4.426763503373093e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 2       2 3   "
                              "-5.352581309875108e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 0    "
                              "1.593864453718162e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 1    "
                              "4.145709149066612e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 2   "
                              "-5.639209206018207e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 3    "
                              "1.576726780865246e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 4   "
                              "-6.823921669935858e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 5    "
                              "2.728030657099649e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 6    "
                              "4.849101303691402e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 7   "
                              "-1.048132260972923e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 8    "
                              "6.745370134377757e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3       2 9   "
                              "-3.536454451630775e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 3      2 10   "
                              "-1.916040882777734e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 0   "
                              "-2.277646296639974e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 1   "
                              "-5.944120868937014e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 2    "
                              "8.084347622615087e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 3   "
                              "-2.258583526759012e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 4    "
                              "9.782633480506241e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 5   "
                              "-3.912513619777572e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 6   "
                              "-6.947830773511784e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 7    "
                              "1.502171795432116e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 8   "
                              "-9.673233375460551e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3       2 9    "
                              "5.068997413252322e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 3      2 10    "
                              "2.746323785359357e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 0   "
                              "-2.065480010246626e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 1   "
                              "-5.176184784101910e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 2    "
                              "6.995720722709600e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 3   "
                              "-1.851981874130438e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 4    "
                              "8.467132979885737e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 5   "
                              "-3.427170472834707e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 6   "
                              "-6.068175346039939e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 7    "
                              "1.296261016708436e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 8   "
                              "-8.386188606718887e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3       2 9    "
                              "4.376268190834145e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 3      2 10    "
                              "2.393366749158289e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 1   "
                              "-1.289598697485621e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 2    "
                              "1.685115800648999e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 3   "
                              "-3.125387959776424e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 4    "
                              "2.041680613739776e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 5   "
                              "-8.815444843999494e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 6   "
                              "-1.527330805382477e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 7    "
                              "3.071950060034302e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 8   "
                              "-2.043744512464694e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 3       2 9    "
                              "1.040082899462996e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 3       2 1    "
                              "4.222429427656757e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 3       2 2   "
                              "-5.391394663410128e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 3       2 4   "
                              "-6.533408977094895e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 3       2 5    "
                              "2.964459475455963e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 3       2 6    "
                              "4.920829119399274e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 3       2 2   "
                              "-5.352581309875108e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 3       2 4   "
                              "-6.486161216250301e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 1    "
                              "5.961941199936058e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 2   "
                              "-6.673225609438826e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 3   "
                              "-6.823921669935858e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 4   "
                              "-1.155886488195298e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 5    "
                              "5.913558666520153e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 6    "
                              "2.301742658267988e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 7    "
                              "2.870199137348181e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 8    "
                              "2.865548953353969e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4       2 9    "
                              "1.077288158496926e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 4      2 10   "
                              "-1.295549003315186e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 1   "
                              "-8.531917225602516e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 2    "
                              "9.567434689579468e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 3    "
                              "9.782633480506240e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 4    "
                              "1.657235531638779e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 5   "
                              "-8.477206560911218e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 6   "
                              "-3.299496388817977e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 7   "
                              "-4.116271372271287e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 8   "
                              "-4.109884916671800e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4       2 9   "
                              "-1.544583845432706e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 4      2 10    "
                              "1.857308648951738e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 1   "
                              "-7.580613009423137e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 2    "
                              "8.328495768830192e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 3    "
                              "8.467132979885735e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 4    "
                              "1.445158440011461e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 5   "
                              "-7.337404555695811e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 6   "
                              "-2.865812163203524e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 7   "
                              "-3.540509740711012e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 8   "
                              "-3.546585198760869e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4       2 9   "
                              "-1.333820104755787e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 4      2 10    "
                              "1.609440180592267e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 1   "
                              "-2.054573374740665e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 2    "
                              "2.070443939231002e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 3    "
                              "2.041680613739776e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 4    "
                              "3.625385721839493e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 5   "
                              "-1.768902579284091e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 6   "
                              "-7.035113844255110e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 7   "
                              "-8.283479859434077e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 8   "
                              "-8.442497964066153e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4       2 9   "
                              "-3.180643948960400e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 4      2 10    "
                              "3.905523455266485e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4       2 2   "
                              "-6.758362232961109e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4       2 3   "
                              "-6.533408977094892e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4       2 4   "
                              "-1.190083045446682e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4       2 5    "
                              "5.651696464201535e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4       2 6    "
                              "2.269296754460614e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4       2 7    "
                              "2.646964005816335e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4       2 8    "
                              "2.721165044127963e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4       2 9    "
                              "1.016320947901377e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 4      2 10   "
                              "-1.255070378117837e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 4       2 3   "
                              "-6.486161216250302e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 4       2 4   "
                              "-1.147919633989248e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 4       2 5    "
                              "5.623506490270192e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 4       2 6    "
                              "2.236000088444864e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5       2 2   "
                              "-2.535838957258822e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5       2 3    "
                              "2.728030657099652e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5       2 4    "
                              "5.913558666520153e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5       2 5    "
                              "5.419833292430022e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5       2 6   "
                              "-7.211522100071873e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5       2 7   "
                              "-1.380397524981230e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5       2 8    "
                              "1.250256146411241e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5       2 9   "
                              "-3.250758670644857e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 5      2 10    "
                              "5.151207797076076e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 1    "
                              "1.211427258319004e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 2    "
                              "3.654065328147384e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 3   "
                              "-3.912513619777576e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 4   "
                              "-8.477206560911219e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 5   "
                              "-7.769467515177771e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 6    "
                              "1.033792629319848e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 7    "
                              "1.978425133722318e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 8   "
                              "-1.788457856768198e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5       2 9    "
                              "4.659399776798742e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 5      2 10   "
                              "-7.386862262533085e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 1    "
                              "1.060793256915104e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 2    "
                              "2.926922722445260e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 3   "
                              "-3.427170472834701e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 4   "
                              "-7.337404555695811e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 5   "
                              "-6.709519553524940e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 6    "
                              "8.954029842646137e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 7    "
                              "1.711338113344460e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 8   "
                              "-1.594447283167990e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5       2 9    "
                              "4.025856583911209e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 5      2 10   "
                              "-6.350926836375669e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 5       2 3   "
                              "-8.815444843999508e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 5       2 4   "
                              "-1.768902579284092e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 5       2 5   "
                              "-1.597955717538821e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 5       2 6    "
                              "2.166537252638215e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 5       2 7    "
                              "4.103989810285592e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 5       2 8   "
                              "-4.382514738828239e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 5       2 9    "
                              "9.607676452642735e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 5      2 10   "
                              "-1.480250734659058e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 5       2 3    "
                              "2.964459475455964e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 5       2 4    "
                              "5.651696464201535e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 5       2 5    "
                              "5.069600097409608e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 5       2 6   "
                              "-6.938272899615851e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 5       2 7   "
                              "-1.296239042406352e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 5       2 8    "
                              "1.417263682290590e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 5       2 9   "
                              "-3.039926879078881e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 5       2 4    "
                              "5.623506490270192e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 5       2 5    "
                              "5.082349703374114e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 5       2 6   "
                              "-6.885951321223108e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 5       2 7   "
                              "-1.309457270819692e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6       2 1    "
                              "1.432939737914284e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6       2 3    "
                              "4.849101303691400e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6       2 4    "
                              "2.301742658267989e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6       2 5   "
                              "-7.211522100071873e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6       2 6   "
                              "-9.718822990700186e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6       2 7    "
                              "3.123201294521006e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6       2 8    "
                              "2.908005904154446e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6       2 9    "
                              "3.631155807796803e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 6      2 10   "
                              "-5.280869161006566e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6       2 1   "
                              "-2.056720199563075e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6       2 3   "
                              "-6.947830773511784e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6       2 4   "
                              "-3.299496388817977e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6       2 5    "
                              "1.033792629319848e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6       2 6    "
                              "1.393234284190626e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6       2 7   "
                              "-4.476962838586570e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6       2 8   "
                              "-4.167697664619744e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6       2 9   "
                              "-5.205931828052437e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 6      2 10    "
                              "7.569898260501175e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6       2 1   "
                              "-1.767087056018783e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6       2 3   "
                              "-6.068175346039928e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6       2 4   "
                              "-2.865812163203524e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6       2 5    "
                              "8.954029842646137e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6       2 6    "
                              "1.208132957006071e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6       2 7   "
                              "-3.875995811715730e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6       2 8   "
                              "-3.624451798408941e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6       2 9   "
                              "-4.493517296103261e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 6      2 10    "
                              "6.571034680643298e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 6       2 3   "
                              "-1.527330805382477e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 6       2 4   "
                              "-7.035113844255110e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 6       2 5    "
                              "2.166537252638214e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 6       2 6    "
                              "2.941473379694497e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 6       2 7   "
                              "-9.352934507908498e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 6       2 8   "
                              "-8.942954506089153e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 6       2 9   "
                              "-1.068556733649410e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 6      2 10    "
                              "1.607878057938480e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 6       2 3    "
                              "4.920829119399272e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 6       2 4    "
                              "2.269296754460616e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 6       2 5   "
                              "-6.938272899615851e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 6       2 6   "
                              "-9.457325127999056e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 6       2 7    "
                              "2.984426762471245e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 6       2 8    "
                              "2.874383580633809e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 6       2 9    "
                              "3.403319184508239e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 6      2 10   "
                              "-5.173506711398747e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 6       2 4    "
                              "2.236000088444864e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 6       2 5   "
                              "-6.885951321223110e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 6       2 6   "
                              "-9.345209158285762e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 6       2 7    "
                              "2.975652051047465e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 6       2 8    "
                              "2.851784702994151e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7       2 2   "
                              "-2.419368543566544e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7       2 3   "
                              "-1.048132260972921e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7       2 4    "
                              "2.870199137348181e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7       2 5   "
                              "-1.380397524981229e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7       2 6    "
                              "3.123201294521006e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7       2 7    "
                              "1.235901989280693e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7       2 8   "
                              "-1.925829476233712e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7       2 9   "
                              "-1.426049970748351e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 7      2 10    "
                              "2.592367564713680e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7       2 2    "
                              "3.469664676211196e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7       2 3    "
                              "1.502171795432107e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7       2 4   "
                              "-4.116271372271284e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7       2 5    "
                              "1.978425133722319e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7       2 6   "
                              "-4.476962838586570e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7       2 7   "
                              "-1.771635016292810e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7       2 8    "
                              "2.760488285775069e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7       2 9    "
                              "2.043755758868814e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 7      2 10   "
                              "-3.716483845585240e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7       2 2    "
                              "2.995335256302941e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7       2 3    "
                              "1.296261016708435e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7       2 4   "
                              "-3.540509740711011e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7       2 5    "
                              "1.711338113344461e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7       2 6   "
                              "-3.875995811715730e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7       2 7   "
                              "-1.532077784662384e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7       2 8    "
                              "2.393694254451009e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7       2 9    "
                              "1.772110467094541e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 7      2 10   "
                              "-3.216256790135352e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 7       2 3    "
                              "3.071950060034298e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 7       2 4   "
                              "-8.283479859434092e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 7       2 5    "
                              "4.103989810285595e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 7       2 6   "
                              "-9.352934507908498e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 7       2 7   "
                              "-3.675033968435040e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 7       2 8    "
                              "5.823096020486472e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 7       2 9    "
                              "4.300596499617232e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 7      2 10   "
                              "-7.750490431838680e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 7       2 4    "
                              "2.646964005816332e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 7       2 5   "
                              "-1.296239042406352e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 7       2 6    "
                              "2.984426762471245e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 7       2 7    "
                              "1.169132270796727e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 7       2 8   "
                              "-1.865051437263648e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 7       2 9   "
                              "-1.363257048405138e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 7      2 10    "
                              "2.480996473328200e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 7       2 5   "
                              "-1.309457270819692e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 7       2 6    "
                              "2.975652051047465e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 7       2 7    "
                              "1.169202560602880e+01 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 7       2 8   "
                              "-1.853176606030544e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 7       2 9   "
                              "-1.372848328346659e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8       2 2   "
                              "-1.860265805809937e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8       2 3    "
                              "6.745370134377756e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8       2 4    "
                              "2.865548953353972e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8       2 5    "
                              "1.250256146411238e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8       2 6    "
                              "2.908005904154445e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8       2 7   "
                              "-1.925829476233712e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8       2 8   "
                              "-1.913770103311484e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8       2 9    "
                              "1.508134375274986e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 8      2 10   "
                              "-2.234716112660251e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8       2 2    "
                              "2.667446337745305e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8       2 3   "
                              "-9.673233375460489e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8       2 4   "
                              "-4.109884916671795e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8       2 5   "
                              "-1.788457856768198e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8       2 6   "
                              "-4.167697664619744e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8       2 7    "
                              "2.760488285775069e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8       2 8    "
                              "2.743215841089073e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8       2 9   "
                              "-2.161757407298414e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 8      2 10    "
                              "3.201621910811511e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8       2 2    "
                              "2.299637577010140e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8       2 3   "
                              "-8.386188606718846e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8       2 4   "
                              "-3.546585198760870e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8       2 5   "
                              "-1.594447283167997e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8       2 6   "
                              "-3.624451798408940e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8       2 7    "
                              "2.393694254451009e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8       2 8    "
                              "2.380772494023351e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8       2 9   "
                              "-1.873508452455056e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 8      2 10    "
                              "2.792924338397143e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 8       2 3   "
                              "-2.043744512464701e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 8       2 4   "
                              "-8.442497964066149e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 8       2 5   "
                              "-4.382514738828238e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 8       2 6   "
                              "-8.942954506089150e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 8       2 7    "
                              "5.823096020486471e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 8       2 8    "
                              "5.818523303454713e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 8       2 9   "
                              "-4.544199297801071e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 8      2 10    "
                              "6.985326764982149e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 8       2 4    "
                              "2.721165044127964e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 8       2 5    "
                              "1.417263682290590e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 8       2 6    "
                              "2.874383580633809e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 8       2 7   "
                              "-1.865051437263647e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 8       2 8   "
                              "-1.869279983043622e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 8       2 9    "
                              "1.452317475658280e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 8      2 10   "
                              "-2.238895790085909e+01 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 8       2 6    "
                              "2.851784702994152e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 8       2 7   "
                              "-1.853176606030544e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 8       2 8   "
                              "-1.851225908332972e+01 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 8       2 9    "
                              "1.446434211493662e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 8      2 10   "
                              "-2.239840669708085e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9       2 2   "
                              "-1.209869728622871e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9       2 3   "
                              "-3.536454451630782e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9       2 4    "
                              "1.077288158496927e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9       2 5   "
                              "-3.250758670644853e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9       2 6    "
                              "3.631155807796803e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9       2 7   "
                              "-1.426049970748351e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9       2 8    "
                              "1.508134375274986e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9       2 9    "
                              "1.559282784835807e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 0       2 9      2 10    "
                              "6.479802931818537e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9       2 2    "
                              "1.734837092442033e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9       2 3    "
                              "5.068997413252322e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9       2 4   "
                              "-1.544583845432706e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9       2 5    "
                              "4.659399776798742e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9       2 6   "
                              "-5.205931828052437e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9       2 7    "
                              "2.043755758868814e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9       2 8   "
                              "-2.161757407298414e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9       2 9   "
                              "-2.235106203535502e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 2       2 9      2 10   "
                              "-9.286708085494013e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9       2 2    "
                              "1.497673462499840e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9       2 3    "
                              "4.376268190834135e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9       2 4   "
                              "-1.333820104755788e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9       2 5    "
                              "4.025856583911219e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9       2 6   "
                              "-4.493517296103261e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9       2 7    "
                              "1.772110467094541e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9       2 8   "
                              "-1.873508452455056e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9       2 9   "
                              "-1.933943941591720e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 4       2 9      2 10   "
                              "-8.068081874153683e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 9       2 3    "
                              "1.040082899462996e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 9       2 4   "
                              "-3.180643948960399e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 9       2 5    "
                              "9.607676452642735e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 9       2 6   "
                              "-1.068556733649410e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 9       2 7    "
                              "4.300596499617232e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 9       2 8   "
                              "-4.544199297801072e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 9       2 9   "
                              "-4.650551535095815e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 6       2 9      2 10   "
                              "-1.979305449693641e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 9       2 4    "
                              "1.016320947901377e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 9       2 5   "
                              "-3.039926879078884e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 9       2 6    "
                              "3.403319184508238e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 9       2 7   "
                              "-1.363257048405137e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 9       2 8    "
                              "1.452317475658280e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 9       2 9    "
                              "1.478825902886419e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 8       2 9      2 10    "
                              "6.329024260738756e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 9       2 7   "
                              "-1.372848328346659e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 9       2 8    "
                              "1.446434211493662e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10       2 9       2 9    "
                              "1.480443351584376e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10       2 2    "
                              "1.230416139854121e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10       2 3   "
                              "-1.916040882777733e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10       2 4   "
                              "-1.295549003315186e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10       2 5    "
                              "5.151207797076060e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10       2 6   "
                              "-5.280869161006563e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10       2 7    "
                              "2.592367564713682e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10       2 8   "
                              "-2.234716112660251e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10       2 9    "
                              "6.479802931818537e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 0      2 10      2 10   "
                              "-2.525337354818048e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10       2 2   "
                              "-1.764152790431516e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10       2 3    "
                              "2.746323785359355e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10       2 4    "
                              "1.857308648951737e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10       2 5   "
                              "-7.386862262533099e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10       2 6    "
                              "7.569898260501174e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10       2 7   "
                              "-3.716483845585240e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10       2 8    "
                              "3.201621910811511e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10       2 9   "
                              "-9.286708085494014e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 2      2 10      2 10    "
                              "3.619338310811442e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10       2 2   "
                              "-1.525628551227811e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10       2 3    "
                              "2.393366749158286e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10       2 4    "
                              "1.609440180592267e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10       2 5   "
                              "-6.350926836375668e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10       2 6    "
                              "6.571034680643298e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10       2 7   "
                              "-3.216256790135352e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10       2 8    "
                              "2.792924338397144e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10       2 9   "
                              "-8.068081874153683e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 4      2 10      2 10    "
                              "3.145579135103967e+03 \n";
    integralFileTwoBodyFAD << "       1 2       1 6      2 10       2 4    "
                              "3.905523455266485e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6      2 10       2 5   "
                              "-1.480250734659062e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6      2 10       2 6    "
                              "1.607878057938481e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6      2 10       2 7   "
                              "-7.750490431838682e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 6      2 10       2 8    "
                              "6.985326764982152e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6      2 10       2 9   "
                              "-1.979305449693642e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 6      2 10      2 10    "
                              "7.737309656852620e+02 \n";
    integralFileTwoBodyFAD << "       1 2       1 8      2 10       2 4   "
                              "-1.255070378117838e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8      2 10       2 6   "
                              "-5.173506711398749e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8      2 10       2 7    "
                              "2.480996473328199e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8      2 10       2 8   "
                              "-2.238895790085907e+01 \n";
    integralFileTwoBodyFAD << "       1 2       1 8      2 10       2 9    "
                              "6.329024260738756e+00 \n";
    integralFileTwoBodyFAD << "       1 2       1 8      2 10      2 10   "
                              "-2.484293963254501e+02 \n";
    integralFileTwoBodyFAD << "       1 2      1 10      2 10       2 8   "
                              "-2.239840669708084e+00 \n";
    integralFileTwoBodyFAD << "       1 2      1 10      2 10      2 10   "
                              "-2.467054026491810e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 0       2 0   "
                              "-2.633211657127596e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 0       2 1    "
                              "3.544633365081027e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 0       2 2   "
                              "-2.625956126994467e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 0       2 3   "
                              "-1.781442696689963e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 0       2 0   "
                              "-5.289281313999776e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 0       2 1    "
                              "6.956509208444363e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 0       2 2   "
                              "-5.136902036282060e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 0       2 3   "
                              "-3.389481107765617e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 0       2 0   "
                              "-3.316175259432389e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 0       2 1    "
                              "4.594621739182922e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 0       2 2   "
                              "-3.417105912761267e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 0       2 3   "
                              "-2.350738657084364e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 0       2 1   "
                              "-1.548304605493568e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 0       2 2    "
                              "1.193959336901428e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 0       2 1   "
                              "-6.555720287020158e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 1       2 0    "
                              "3.544633365081026e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 1       2 2   "
                              "-5.013246559842776e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 1       2 3   "
                              "-4.473105088136509e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 1       2 4   "
                              "-6.547100969860336e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 1       2 6   "
                              "-1.528354157961985e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 1       2 0    "
                              "6.956509208444363e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 1       2 2   "
                              "-9.837206083207317e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 1       2 3   "
                              "-8.749565696460725e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 1       2 4   "
                              "-1.262490162877374e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 1       2 5    "
                              "1.783557507772168e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 1       2 6   "
                              "-3.016501186335827e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 1       2 7    "
                              "1.054647329643900e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 1       2 8   "
                              "-1.337218133533638e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 1       2 0    "
                              "4.594621739182922e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 1       2 2   "
                              "-6.498731144299226e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 1       2 3   "
                              "-5.820462306712188e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 1       2 4   "
                              "-8.574476810455801e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 1       2 5    "
                              "1.193383332023026e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 1       2 6   "
                              "-1.978363485299688e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 1       2 0   "
                              "-1.548304605493568e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 1       2 1   "
                              "-1.052824373302287e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 1       2 2    "
                              "2.191766868201736e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 1       2 3    "
                              "2.033046895459614e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 1       2 4    "
                              "3.209946277326264e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 1       2 0   "
                              "-6.555720287020159e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 1       2 2    "
                              "9.280444835313176e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 0   "
                              "-2.625956126994464e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 1   "
                              "-5.013246559842777e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 2    "
                              "3.486206322710816e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 3    "
                              "6.053411383368234e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 4    "
                              "7.198349927955805e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 5    "
                              "2.546156025363763e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 7    "
                              "2.591506863007378e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 8    "
                              "1.990500127001242e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2       2 9    "
                              "1.296047570331961e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 2      2 10   "
                              "-1.320012570275534e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 0   "
                              "-5.136902036282060e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 1   "
                              "-9.837206083207318e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 2    "
                              "6.256910385125215e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 3    "
                              "1.187631097640627e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 4    "
                              "1.408106185045677e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 5    "
                              "5.260471825681455e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 7    "
                              "5.094336321774373e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 8    "
                              "3.914594051023573e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2       2 9    "
                              "2.546996591532522e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 2      2 10   "
                              "-2.591293632279440e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 0   "
                              "-3.417105912761269e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 1   "
                              "-6.498731144299226e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 2    "
                              "4.978911547778019e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 3    "
                              "7.847674701878742e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 4    "
                              "9.363634812346245e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 5    "
                              "3.200745813643807e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 7    "
                              "3.358013897913122e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 8    "
                              "2.576579352025349e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2       2 9    "
                              "1.678859903465051e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 2      2 10   "
                              "-1.711175652333371e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 2       2 0    "
                              "1.193959336901429e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 2       2 1    "
                              "2.191766868201736e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 2       2 2   "
                              "-3.149375302594244e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 2       2 3   "
                              "-2.648950877321754e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 2       2 4   "
                              "-3.262191804448531e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 2       2 7   "
                              "-1.125542857857236e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 2       2 1    "
                              "9.280444835313176e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 2       2 3   "
                              "-1.121720773849247e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 2       2 4   "
                              "-1.340400183023993e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 0   "
                              "-1.781442696689939e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 1   "
                              "-4.473105088136506e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 2    "
                              "6.053411383368234e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 3   "
                              "-1.620015215629238e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 4    "
                              "7.326486073470033e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 5   "
                              "-2.957236128215965e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 6   "
                              "-5.246678329675058e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 7    "
                              "1.122598623169609e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 8   "
                              "-7.251075798011829e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3       2 9    "
                              "3.789217507298867e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 3      2 10    "
                              "2.069022110495707e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 0   "
                              "-3.389481107765642e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 1   "
                              "-8.749565696460731e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 2    "
                              "1.187631097640626e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 3   "
                              "-3.263852339178073e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 4    "
                              "1.437203142966596e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 5   "
                              "-5.770376969151326e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 6   "
                              "-1.023288253418356e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 7    "
                              "2.204655561041955e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 8   "
                              "-1.422041872431124e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3       2 9    "
                              "7.440827739506768e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 3      2 10    "
                              "4.042822084172319e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 0   "
                              "-2.350738657084360e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 1   "
                              "-5.820462306712191e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 2    "
                              "7.847674701878742e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 3   "
                              "-2.034217450612956e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 4    "
                              "9.498941090852404e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 5   "
                              "-3.863105522499197e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 6   "
                              "-6.827676857536898e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 7    "
                              "1.452537613796562e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 8   "
                              "-9.415227099173276e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3       2 9    "
                              "4.904632855030777e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 3      2 10    "
                              "2.691417579637127e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 1    "
                              "2.033046895459614e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 2   "
                              "-2.648950877321754e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 3    "
                              "4.752201350175793e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 4   "
                              "-3.209341015703793e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 5    "
                              "1.395360089645934e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 6    "
                              "2.396859159788275e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 7   "
                              "-4.816566679604802e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 8    "
                              "3.221609331654029e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 3       2 9   "
                              "-1.631812734709111e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 3       2 2   "
                              "-1.121720773849247e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 3       2 4   "
                              "-1.358501369865259e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 3       2 6    "
                              "1.000129368049811e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 1   "
                              "-6.547100969860321e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 2    "
                              "7.198349927955805e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 3    "
                              "7.326486073470032e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 4    "
                              "1.248631492851935e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 5   "
                              "-6.349332585316283e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 6   "
                              "-2.478479621974622e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 7   "
                              "-3.064871319823030e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 8   "
                              "-3.068365192265271e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4       2 9   "
                              "-1.154350779134184e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 4      2 10    "
                              "1.392298577704566e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 1   "
                              "-1.262490162877369e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 2    "
                              "1.408106185045676e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 3    "
                              "1.437203142966596e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 4    "
                              "2.440401782737295e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 5   "
                              "-1.245400522585288e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 6   "
                              "-4.852426297375032e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 7   "
                              "-6.037180588841880e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 8   "
                              "-6.033765429213899e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4       2 9   "
                              "-2.267781613374260e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 4      2 10    "
                              "2.729641822193031e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 1   "
                              "-8.574476810455772e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 2    "
                              "9.363634812346251e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 3    "
                              "9.498941090852404e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 4    "
                              "1.625842713109470e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 5   "
                              "-8.231351574909822e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 6   "
                              "-3.219045746185515e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 7   "
                              "-3.964321159021027e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 8   "
                              "-3.975455982742326e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4       2 9   "
                              "-1.495263735738211e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 4      2 10    "
                              "1.806368196341459e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 1    "
                              "3.209946277326268e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 2   "
                              "-3.262191804448532e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 3   "
                              "-3.209341015703792e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 4   "
                              "-5.715926739403926e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 5    "
                              "2.779598137705229e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 6    "
                              "1.106400789464097e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 7    "
                              "1.304669400881714e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 8    "
                              "1.330679042005104e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4       2 9    "
                              "5.002264664759646e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 4      2 10   "
                              "-6.142198838335274e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 4       2 2   "
                              "-1.340400183023993e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 4       2 3   "
                              "-1.358501369865259e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 4       2 4   "
                              "-2.328385456519784e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 4       2 5    "
                              "1.178916342433936e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 4       2 6    "
                              "4.625288618906213e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5       2 2    "
                              "2.546156025363763e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5       2 3   "
                              "-2.957236128215959e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5       2 4   "
                              "-6.349332585316283e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5       2 5   "
                              "-5.808308224762370e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5       2 6    "
                              "7.747258227388839e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5       2 7    "
                              "1.481572049056897e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5       2 8   "
                              "-1.376902317603778e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5       2 9    "
                              "3.485434155131458e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 5      2 10   "
                              "-5.499337648114057e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 1    "
                              "1.783557507772182e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 2    "
                              "5.260471825681455e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 3   "
                              "-5.770376969151337e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 4   "
                              "-1.245400522585288e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 5   "
                              "-1.140635730787463e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 6    "
                              "1.519081249011605e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 7    "
                              "2.905649495592583e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 8   "
                              "-2.648832036383267e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5       2 9    "
                              "6.841163882065715e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 5      2 10   "
                              "-1.083186074253887e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 1    "
                              "1.193383332022997e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 2    "
                              "3.200745813643778e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 3   "
                              "-3.863105522499191e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 4   "
                              "-8.231351574909822e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 5   "
                              "-7.520633571223505e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 6    "
                              "1.004750511298486e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 7    "
                              "1.918985185655001e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 8   "
                              "-1.805567291700102e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5       2 9    "
                              "4.513174695465789e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 5      2 10   "
                              "-7.108706853343778e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 5       2 3    "
                              "1.395360089645934e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 5       2 4    "
                              "2.779598137705229e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 5       2 5    "
                              "2.509255609007438e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 5       2 6   "
                              "-3.405310144163746e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 5       2 7   "
                              "-6.434589495140034e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 5       2 8    "
                              "6.840243717938977e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 5       2 9   "
                              "-1.507482553964327e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 5      2 10    "
                              "2.327909692629900e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 5       2 4    "
                              "1.178916342433936e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 5       2 5    "
                              "1.075395185438978e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 5       2 6   "
                              "-1.439389276842474e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 5       2 7   "
                              "-2.768707399622573e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6       2 1   "
                              "-1.528354157961956e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6       2 3   "
                              "-5.246678329675058e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6       2 4   "
                              "-2.478479621974622e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6       2 5    "
                              "7.747258227388841e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6       2 6    "
                              "1.045076850443366e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6       2 7   "
                              "-3.354159511611819e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6       2 8   "
                              "-3.134849417383397e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6       2 9   "
                              "-3.889522473432193e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 6      2 10    "
                              "5.683666399448375e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6       2 1   "
                              "-3.016501186335784e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6       2 3   "
                              "-1.023288253418356e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6       2 4   "
                              "-4.852426297375032e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6       2 5    "
                              "1.519081249011605e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6       2 6    "
                              "2.047991212864563e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6       2 7   "
                              "-6.577522044927815e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6       2 8   "
                              "-6.131013231942054e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6       2 9   "
                              "-7.642182193507307e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 6      2 10    "
                              "1.113061370597894e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6       2 1   "
                              "-1.978363485299688e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6       2 3   "
                              "-6.827676857536899e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6       2 4   "
                              "-3.219045746185515e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6       2 5    "
                              "1.004750511298486e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6       2 6    "
                              "1.356266944253072e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6       2 7   "
                              "-4.348436915246994e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6       2 8   "
                              "-4.072545088132770e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6       2 9   "
                              "-5.036452631835039e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 6      2 10    "
                              "7.379221537115292e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 6       2 3    "
                              "2.396859159788276e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 6       2 4    "
                              "1.106400789464097e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 6       2 5   "
                              "-3.405310144163745e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 6       2 6   "
                              "-4.625384216280532e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 6       2 7    "
                              "1.469152080985388e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 6       2 8    "
                              "1.404713045722655e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 6       2 9    "
                              "1.679625509592283e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 6      2 10   "
                              "-2.527832710654951e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 6       2 3    "
                              "1.000129368049812e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 6       2 4    "
                              "4.625288618906215e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 6       2 5   "
                              "-1.439389276842474e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 6       2 6   "
                              "-1.943779318542299e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 6       2 7    "
                              "6.239953534808506e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 6       2 8    "
                              "5.897004744747056e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 6      2 10   "
                              "-1.060383707821506e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7       2 2    "
                              "2.591506863007378e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7       2 3    "
                              "1.122598623169610e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7       2 4   "
                              "-3.064871319823030e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7       2 5    "
                              "1.481572049056897e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7       2 6   "
                              "-3.354159511611819e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7       2 7   "
                              "-1.326046253590376e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7       2 8    "
                              "2.070948914317344e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7       2 9    "
                              "1.533730031004695e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 7      2 10   "
                              "-2.782946462827163e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 1    "
                              "1.054647329643900e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 2    "
                              "5.094336321774359e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 3    "
                              "2.204655561041960e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 4   "
                              "-6.037180588841891e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 5    "
                              "2.905649495592580e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 6   "
                              "-6.577522044927816e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 7   "
                              "-2.601991450290603e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 8    "
                              "4.057574379037761e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7       2 9    "
                              "3.003672196485837e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 7      2 10   "
                              "-5.459848707939211e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7       2 2    "
                              "3.358013897913129e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7       2 3    "
                              "1.452537613796560e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7       2 4   "
                              "-3.964321159021029e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7       2 5    "
                              "1.918985185655002e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7       2 6   "
                              "-4.348436915246994e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7       2 7   "
                              "-1.718113719975473e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7       2 8    "
                              "2.686979812925845e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7       2 9    "
                              "1.988655440864239e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 7      2 10   "
                              "-3.608016257964326e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7       2 2   "
                              "-1.125542857857247e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7       2 3   "
                              "-4.816566679604802e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7       2 4    "
                              "1.304669400881716e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7       2 5   "
                              "-6.434589495140037e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7       2 6    "
                              "1.469152080985388e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7       2 7    "
                              "5.771349794646538e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7       2 8   "
                              "-9.149074069379840e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7       2 9   "
                              "-6.742517741361073e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 7      2 10    "
                              "1.218472311167989e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 7       2 5   "
                              "-2.768707399622573e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 7       2 6    "
                              "6.239953534808506e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 7       2 7    "
                              "2.462404879840945e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 7       2 8   "
                              "-3.864053571595488e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 7       2 9   "
                              "-2.882133528804209e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8       2 2    "
                              "1.990500127001243e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8       2 3   "
                              "-7.251075798011847e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8       2 4   "
                              "-3.068365192265273e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8       2 5   "
                              "-1.376902317603778e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8       2 6   "
                              "-3.134849417383397e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8       2 7    "
                              "2.070948914317344e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8       2 8    "
                              "2.059427107508353e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8       2 9   "
                              "-1.621077623892800e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 8      2 10    "
                              "2.415562839422610e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 1   "
                              "-1.337218133533625e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 2    "
                              "3.914594051023573e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 3   "
                              "-1.422041872431121e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 4   "
                              "-6.033765429213898e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 5   "
                              "-2.648832036383267e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 6   "
                              "-6.131013231942054e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 7    "
                              "4.057574379037761e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 8    "
                              "4.033268098305828e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8       2 9   "
                              "-3.176976432252673e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 8      2 10    "
                              "4.713539105109805e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8       2 2    "
                              "2.576579352025334e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8       2 3   "
                              "-9.415227099173276e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8       2 4   "
                              "-3.975455982742329e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8       2 5   "
                              "-1.805567291700101e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8       2 6   "
                              "-4.072545088132769e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8       2 7    "
                              "2.686979812925845e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8       2 8    "
                              "2.673363584959902e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8       2 9   "
                              "-2.102601543823059e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 8      2 10    "
                              "3.141131878579669e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 8       2 3    "
                              "3.221609331654015e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 8       2 4    "
                              "1.330679042005102e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 8       2 5    "
                              "6.840243717938991e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 8       2 6    "
                              "1.404713045722655e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 8       2 7   "
                              "-9.149074069379840e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 8       2 8   "
                              "-9.145112929589288e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 8       2 9    "
                              "7.137817176133412e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 8      2 10   "
                              "-1.095365014311528e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 8       2 6    "
                              "5.897004744747056e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 8       2 7   "
                              "-3.864053571595488e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 8       2 8   "
                              "-3.845511104588638e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 8       2 9    "
                              "3.023552700596141e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 8      2 10   "
                              "-4.610521400014510e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9       2 2    "
                              "1.296047570331961e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9       2 3    "
                              "3.789217507298869e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9       2 4   "
                              "-1.154350779134184e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9       2 5    "
                              "3.485434155131459e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9       2 6   "
                              "-3.889522473432192e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9       2 7    "
                              "1.533730031004695e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9       2 8   "
                              "-1.621077623892800e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9       2 9   "
                              "-1.673834988055630e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 1       2 9      2 10   "
                              "-6.979851870801234e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9       2 2    "
                              "2.546996591532507e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9       2 3    "
                              "7.440827739506768e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9       2 4   "
                              "-2.267781613374263e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9       2 5    "
                              "6.841163882065720e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9       2 6   "
                              "-7.642182193507313e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9       2 7    "
                              "3.003672196485835e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9       2 8   "
                              "-3.176976432252673e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9       2 9   "
                              "-3.283157948708549e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3       2 9      2 10   "
                              "-1.365700039419134e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9       2 2    "
                              "1.678859903465043e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9       2 3    "
                              "4.904632855030791e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9       2 4   "
                              "-1.495263735738213e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9       2 5    "
                              "4.513174695465784e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9       2 6   "
                              "-5.036452631835039e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9       2 7    "
                              "1.988655440864239e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9       2 8   "
                              "-2.102601543823059e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9       2 9   "
                              "-2.169100351846719e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 5       2 9      2 10   "
                              "-9.061377781098233e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 9       2 3   "
                              "-1.631812734709108e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 9       2 4    "
                              "5.002264664759639e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 9       2 5   "
                              "-1.507482553964329e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 9       2 6    "
                              "1.679625509592283e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 9       2 7   "
                              "-6.742517741361075e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 9       2 8    "
                              "7.137817176133413e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 9       2 9    "
                              "7.301264573074578e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7       2 9      2 10    "
                              "3.106728481330241e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 9       2 7   "
                              "-2.882133528804209e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 9       2 8    "
                              "3.023552700596141e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 9       2 9    "
                              "3.115214987854690e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 9       2 9      2 10    "
                              "1.313011779560913e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10       2 2   "
                              "-1.320012570275535e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10       2 3    "
                              "2.069022110495703e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10       2 4    "
                              "1.392298577704565e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10       2 5   "
                              "-5.499337648114057e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10       2 6    "
                              "5.683666399448372e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10       2 7   "
                              "-2.782946462827163e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10       2 8    "
                              "2.415562839422610e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10       2 9   "
                              "-6.979851870801234e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 1      2 10      2 10    "
                              "2.720876018621629e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10       2 2   "
                              "-2.591293632279444e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10       2 3    "
                              "4.042822084172317e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10       2 4    "
                              "2.729641822193033e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10       2 5   "
                              "-1.083186074253888e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10       2 6    "
                              "1.113061370597894e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10       2 7   "
                              "-5.459848707939208e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10       2 8    "
                              "4.713539105109807e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10       2 9   "
                              "-1.365700039419134e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 3      2 10      2 10    "
                              "5.323361576220868e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10       2 2   "
                              "-1.711175652333362e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10       2 3    "
                              "2.691417579637133e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10       2 4    "
                              "1.806368196341460e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10       2 5   "
                              "-7.108706853343776e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10       2 6    "
                              "7.379221537115291e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10       2 7   "
                              "-3.608016257964326e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10       2 8    "
                              "3.141131878579668e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10       2 9   "
                              "-9.061377781098231e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 5      2 10      2 10    "
                              "3.533710787834433e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 7      2 10       2 4   "
                              "-6.142198838335274e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7      2 10       2 5    "
                              "2.327909692629905e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 7      2 10       2 6   "
                              "-2.527832710654953e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7      2 10       2 7    "
                              "1.218472311167990e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7      2 10       2 8   "
                              "-1.095365014311529e+02 \n";
    integralFileTwoBodyFAD << "       1 3       1 7      2 10       2 9    "
                              "3.106728481330242e+01 \n";
    integralFileTwoBodyFAD << "       1 3       1 7      2 10      2 10   "
                              "-1.215325755672726e+03 \n";
    integralFileTwoBodyFAD << "       1 3       1 9      2 10       2 6   "
                              "-1.060383707821505e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9      2 10       2 8   "
                              "-4.610521400014510e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9      2 10       2 9    "
                              "1.313011779560913e+00 \n";
    integralFileTwoBodyFAD << "       1 3       1 9      2 10      2 10   "
                              "-5.111371542171214e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 0       2 1   "
                              "-1.392468260880556e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 0       2 2    "
                              "1.046968311155214e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 0       2 0   "
                              "-3.012580490735589e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 0       2 1    "
                              "4.096365102835440e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 0       2 2   "
                              "-3.038853937095197e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 0       2 3   "
                              "-2.065480010246668e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 0       2 0   "
                              "-5.784781254106729e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 0       2 1    "
                              "7.708639329175583e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 0       2 2   "
                              "-5.702778323602459e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 0       2 3   "
                              "-3.797077013488432e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 0       2 0   "
                              "-3.372303250186400e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 0       2 1    "
                              "4.901474884730100e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 0       2 2   "
                              "-3.667946726863763e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 0       2 3   "
                              "-2.575222147295247e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 0       2 1   "
                              "-1.330708936734694e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 0       2 2    "
                              "1.042637926696740e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 0       2 1   "
                              "-5.718389350457167e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 1       2 0   "
                              "-1.392468260880556e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 1       2 2    "
                              "1.970386576852509e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 1       2 3    "
                              "1.783532735367434e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 1       2 4    "
                              "2.728398926900654e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 1       2 0    "
                              "4.096365102835439e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 1       2 2   "
                              "-5.793608983847445e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 1       2 3   "
                              "-5.176184784101913e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 1       2 4   "
                              "-7.580613009423192e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 1       2 5    "
                              "1.060793256915089e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 1       2 6   "
                              "-1.767087056018877e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 1       2 0    "
                              "7.708639329175583e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 1       2 2   "
                              "-1.090130448612136e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 1       2 3   "
                              "-9.713386042526157e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 1       2 4   "
                              "-1.407668727732162e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 1       2 5    "
                              "1.980589876100330e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 1       2 6   "
                              "-3.338264160936500e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 1       2 7    "
                              "1.167334264612688e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 1       2 8   "
                              "-1.482398297344683e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 1       2 0    "
                              "4.901474884730099e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 1       2 2   "
                              "-6.933484505246231e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 1       2 3   "
                              "-6.247104701786767e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 1       2 4   "
                              "-9.289958677015017e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 1       2 5    "
                              "1.272510768566127e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 1       2 6   "
                              "-2.106840597584691e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 1       2 0   "
                              "-1.330708936734694e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 1       2 1   "
                              "-1.259122793097369e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 1       2 2    "
                              "1.884671390058140e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 1       2 3    "
                              "1.775379814730941e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 1       2 4    "
                              "2.910704981605819e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 1       2 0   "
                              "-5.718389350457166e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 1       2 2    "
                              "8.094570959130641e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 2       2 0    "
                              "1.046968311155213e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 2       2 1    "
                              "1.970386576852509e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 2       2 2   "
                              "-1.905519886316620e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 2       2 3   "
                              "-2.380460784189811e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 2       2 4   "
                              "-2.868054581076223e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 2       2 7   "
                              "-1.013408184746762e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 0   "
                              "-3.038853937095197e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 1   "
                              "-5.793608983847445e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 2    "
                              "4.171563983873238e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 3    "
                              "6.995720722709600e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 4    "
                              "8.328495768830192e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 5    "
                              "2.926922722445296e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 7    "
                              "2.995335256302905e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 8    "
                              "2.299637577010125e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2       2 9    "
                              "1.497673462499825e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 2      2 10   "
                              "-1.525628551227811e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 0   "
                              "-5.702778323602459e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 1   "
                              "-1.090130448612136e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 2    "
                              "7.297800160218027e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 3    "
                              "1.316161600882796e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 4    "
                              "1.563032940341473e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 5    "
                              "5.729865730116564e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 7    "
                              "5.643149123383671e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 8    "
                              "4.334481872461538e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2       2 9    "
                              "2.821190975353558e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 2      2 10   "
                              "-2.871439969097596e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 0   "
                              "-3.667946726863763e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 1   "
                              "-6.933484505246231e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 2    "
                              "6.093738515021649e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 3    "
                              "8.373557364342490e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 4    "
                              "1.004477670517755e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 5    "
                              "3.253023759749102e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 7    "
                              "3.580764339138660e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 8    "
                              "2.742839545884803e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2       2 9    "
                              "1.789233425378182e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 2      2 10   "
                              "-1.825769883946858e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 2       2 0    "
                              "1.042637926696740e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 2       2 1    "
                              "1.884671390058140e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 2       2 2   "
                              "-3.279620017120488e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 2       2 3   "
                              "-2.278966218560902e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 2       2 4   "
                              "-2.846256946553360e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 2       2 1    "
                              "8.094570959130641e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 2       2 3   "
                              "-9.783918468214843e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 2       2 4   "
                              "-1.120380157150110e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 1    "
                              "1.783532735367432e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 2   "
                              "-2.380460784189811e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 3    "
                              "5.595564612473130e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 4   "
                              "-2.882497074847803e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 5    "
                              "1.193925048101817e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 6    "
                              "2.106563414627799e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 7   "
                              "-4.387073529628132e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 8    "
                              "2.862872488740670e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 3       2 9   "
                              "-1.482495648159447e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 0   "
                              "-2.065480010246647e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 1   "
                              "-5.176184784101913e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 2    "
                              "6.995720722709600e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 3   "
                              "-1.851981874130432e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 4    "
                              "8.467132979885737e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 5   "
                              "-3.427170472834708e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 6   "
                              "-6.068175346039932e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 7    "
                              "1.296261016708435e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 8   "
                              "-8.386188606718845e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3       2 9    "
                              "4.376268190834156e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 3      2 10    "
                              "2.393366749158285e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 0   "
                              "-3.797077013488407e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 1   "
                              "-9.713386042526180e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 2    "
                              "1.316161600882796e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 3   "
                              "-3.564475650432418e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 4    "
                              "1.592823569372923e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 5   "
                              "-6.417432292092322e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 6   "
                              "-1.136513082219092e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 7    "
                              "2.441312809244151e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 8   "
                              "-1.576891164176069e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3       2 9    "
                              "8.240513849039976e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 3      2 10    "
                              "4.488331227771418e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 0   "
                              "-2.575222147295244e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 1   "
                              "-6.247104701786767e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 2    "
                              "8.373557364342490e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 3   "
                              "-2.058523508936851e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 4    "
                              "1.013686331149331e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 5   "
                              "-4.172003642850807e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 6   "
                              "-7.327405851298596e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 7    "
                              "1.544951744941334e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 8   "
                              "-1.007254829333756e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3       2 9    "
                              "5.219917234123769e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 3      2 10    "
                              "2.886987091531486e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 1    "
                              "1.775379814730943e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 2   "
                              "-2.278966218560902e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 3    "
                              "3.263834189772104e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 4   "
                              "-2.762467533221910e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 5    "
                              "1.234450371072650e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 6    "
                              "2.104573333345826e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 7   "
                              "-4.113180474383513e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 8    "
                              "2.785479656651576e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 3       2 9   "
                              "-1.395412567391568e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 3       2 2   "
                              "-9.783918468214845e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 3       2 4   "
                              "-1.184225594763195e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 1    "
                              "2.728398926900643e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 2   "
                              "-2.868054581076223e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 3   "
                              "-2.882497074847803e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 4   "
                              "-4.993945515708541e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 5    "
                              "2.498254853469479e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 6    "
                              "9.827606592886532e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 7    "
                              "1.188386623087853e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 8    "
                              "1.198793496548524e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4       2 9    "
                              "4.517805884059684e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 4      2 10   "
                              "-5.492070086242435e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 1   "
                              "-7.580613009423165e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 2    "
                              "8.328495768830192e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 3    "
                              "8.467132979885737e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 4    "
                              "1.445158440011461e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 5   "
                              "-7.337404555695811e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 6   "
                              "-2.865812163203524e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 7   "
                              "-3.540509740711011e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 8   "
                              "-3.546585198760869e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4       2 9   "
                              "-1.333820104755788e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 4      2 10    "
                              "1.609440180592268e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 1   "
                              "-1.407668727732162e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 2    "
                              "1.563032940341472e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 3    "
                              "1.592823569372923e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 4    "
                              "2.710207957780472e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 5   "
                              "-1.380228532046253e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 6   "
                              "-5.382702020571792e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 7   "
                              "-6.681638221607770e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 8   "
                              "-6.683123660819747e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4       2 9   "
                              "-2.512010284802133e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 4      2 10    "
                              "3.026180852572170e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 1   "
                              "-9.289958677014988e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 2    "
                              "1.004477670517754e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 3    "
                              "1.013686331149331e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 4    "
                              "1.746846391965047e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 5   "
                              "-8.782855812483540e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 6   "
                              "-3.444658859938158e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 7   "
                              "-4.215314041503879e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 8   "
                              "-4.238545057657451e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4       2 9   "
                              "-1.593432574526107e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 4      2 10    "
                              "1.929757887422542e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 1    "
                              "2.910704981605826e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 2   "
                              "-2.846256946553358e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 3   "
                              "-2.762467533221910e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 4   "
                              "-5.006662557819734e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 5    "
                              "2.392495473378873e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 6    "
                              "9.601837804309110e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 7    "
                              "1.106297935906884e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 8    "
                              "1.137894864898774e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4       2 9    "
                              "4.282577880577520e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 4      2 10   "
                              "-5.302150443587372e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 4       2 2   "
                              "-1.120380157150111e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 4       2 3   "
                              "-1.184225594763195e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 4       2 4   "
                              "-1.921375436823272e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 4       2 5    "
                              "1.030261558620388e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 4       2 6    "
                              "3.960782823022999e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 5       2 3    "
                              "1.193925048101817e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 5       2 4    "
                              "2.498254853469479e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 5       2 5    "
                              "2.273767402683973e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 5       2 6   "
                              "-3.052896196016509e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 5       2 7   "
                              "-5.822809404803633e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 5       2 8    "
                              "5.777886299654134e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 5       2 9   "
                              "-1.366111317915933e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 5      2 10    "
                              "2.130802009503947e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 1    "
                              "1.060793256915091e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 2    "
                              "2.926922722445324e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 3   "
                              "-3.427170472834708e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 4   "
                              "-7.337404555695811e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 5   "
                              "-6.709519553524940e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 6    "
                              "8.954029842646137e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 7    "
                              "1.711338113344459e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 8   "
                              "-1.594447283167990e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5       2 9    "
                              "4.025856583911201e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 5      2 10   "
                              "-6.350926836375668e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 1    "
                              "1.980589876100359e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 2    "
                              "5.729865730116678e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 3   "
                              "-6.417432292092344e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 4   "
                              "-1.380228532046253e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 5   "
                              "-1.263353034342018e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 6    "
                              "1.683849949593440e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 7    "
                              "3.219177986311282e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 8   "
                              "-2.955970298514095e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5       2 9    "
                              "7.577926689386558e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 5      2 10   "
                              "-1.198518721277427e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 1    "
                              "1.272510768566127e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 2    "
                              "3.253023759749045e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 3   "
                              "-4.172003642850819e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 4   "
                              "-8.782855812483540e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 5   "
                              "-8.008799078254622e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 6    "
                              "1.072723553863097e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 7    "
                              "2.044489550854508e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 8   "
                              "-1.959727484419260e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5       2 9    "
                              "4.806215534947030e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 5      2 10   "
                              "-7.551003996121437e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 5       2 3    "
                              "1.234450371072650e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 5       2 4    "
                              "2.392495473378872e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 5       2 5    "
                              "2.147640920859369e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 5       2 6   "
                              "-2.935943264637376e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 5       2 7   "
                              "-5.527434010713306e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 5       2 8    "
                              "6.237387858775494e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 5       2 9   "
                              "-1.291638838655702e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 5      2 10    "
                              "1.971049907245196e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 5       2 4    "
                              "1.030261558620388e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 5       2 5    "
                              "9.531551167504967e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 5       2 6   "
                              "-1.252039584713831e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 5       2 7   "
                              "-2.464810671566449e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 6       2 3    "
                              "2.106563414627798e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 6       2 4    "
                              "9.827606592886529e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 6       2 5   "
                              "-3.052896196016509e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 6       2 6   "
                              "-4.128850191326191e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 6       2 7    "
                              "1.320567810911419e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 6       2 8    "
                              "1.246636358465198e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 6       2 9    "
                              "1.520722898415504e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 6      2 10   "
                              "-2.250794069603903e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6       2 1   "
                              "-1.767087056018841e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6       2 3   "
                              "-6.068175346039921e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6       2 4   "
                              "-2.865812163203524e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6       2 5    "
                              "8.954029842646137e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6       2 6    "
                              "1.208132957006071e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6       2 7   "
                              "-3.875995811715730e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6       2 8   "
                              "-3.624451798408941e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6       2 9   "
                              "-4.493517296103259e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 6      2 10    "
                              "6.571034680643295e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6       2 1   "
                              "-3.338264160936557e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6       2 3   "
                              "-1.136513082219093e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6       2 4   "
                              "-5.382702020571794e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6       2 5    "
                              "1.683849949593440e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6       2 6    "
                              "2.270851960051119e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6       2 7   "
                              "-7.289864127542392e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6       2 8   "
                              "-6.802638916863284e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6       2 9   "
                              "-8.464039591527386e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 6      2 10    "
                              "1.234485106730705e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6       2 1   "
                              "-2.106840597584676e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6       2 3   "
                              "-7.327405851298607e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6       2 4   "
                              "-3.444658859938158e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6       2 5    "
                              "1.072723553863097e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6       2 6    "
                              "1.449533889084441e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6       2 7   "
                              "-4.639822205897415e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6       2 8   "
                              "-4.359343171640757e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6       2 9   "
                              "-5.364012338293058e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 6      2 10    "
                              "7.891806030829299e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 6       2 3    "
                              "2.104573333345826e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 6       2 4    "
                              "9.601837804309110e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 6       2 5   "
                              "-2.935943264637376e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 6       2 6   "
                              "-3.999055142196680e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 6       2 7    "
                              "1.265206707577585e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 6       2 8    "
                              "1.222199276527281e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 6       2 9    "
                              "1.436182314881395e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 6      2 10   "
                              "-2.190682304157161e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 6       2 4    "
                              "3.960782823022997e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 6       2 5   "
                              "-1.252039584713831e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 6       2 6   "
                              "-1.677232895230527e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 6       2 7    "
                              "5.462478120731875e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 6       2 8    "
                              "5.074230704757887e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7       2 2   "
                              "-1.013408184746753e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7       2 3   "
                              "-4.387073529628132e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7       2 4    "
                              "1.188386623087853e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7       2 5   "
                              "-5.822809404803633e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7       2 6    "
                              "1.320567810911419e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7       2 7    "
                              "5.207536126860405e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7       2 8   "
                              "-8.182168002205954e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7       2 9   "
                              "-6.061114169503511e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 7      2 10    "
                              "1.094561935310272e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7       2 2    "
                              "2.995335256302941e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7       2 3    "
                              "1.296261016708435e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7       2 4   "
                              "-3.540509740711011e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7       2 5    "
                              "1.711338113344461e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7       2 6   "
                              "-3.875995811715730e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7       2 7   "
                              "-1.532077784662384e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7       2 8    "
                              "2.393694254451009e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7       2 9    "
                              "1.772110467094541e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 7      2 10   "
                              "-3.216256790135353e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 1    "
                              "1.167334264612691e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 2    "
                              "5.643149123383627e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 3    "
                              "2.441312809244151e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 4   "
                              "-6.681638221607774e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 5    "
                              "3.219177986311281e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 6   "
                              "-7.289864127542392e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 7   "
                              "-2.882927797406454e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 8    "
                              "4.498850162058931e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7       2 9    "
                              "3.329629255617322e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 7      2 10   "
                              "-6.050833727024074e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7       2 2    "
                              "3.580764339138653e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7       2 3    "
                              "1.544951744941334e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7       2 4   "
                              "-4.215314041503885e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7       2 5    "
                              "2.044489550854508e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7       2 6   "
                              "-4.639822205897415e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7       2 7   "
                              "-1.831543874908211e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7       2 8    "
                              "2.870595717144795e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7       2 9    "
                              "2.122274753797954e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 7      2 10   "
                              "-3.850153106984001e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 7       2 3   "
                              "-4.113180474383506e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 7       2 4    "
                              "1.106297935906884e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 7       2 5   "
                              "-5.527434010713306e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 7       2 6    "
                              "1.265206707577585e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 7       2 7    "
                              "4.956485384580501e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 7       2 8   "
                              "-7.908438434585075e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 7       2 9   "
                              "-5.825026979372001e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 7      2 10    "
                              "1.048492818872038e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 7       2 5   "
                              "-2.464810671566449e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 7       2 6    "
                              "5.462478120731876e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 7       2 7    "
                              "2.169086033921918e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 7       2 8   "
                              "-3.355787621948815e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 7       2 9   "
                              "-2.542766261503350e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 8       2 3    "
                              "2.862872488740655e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 8       2 4    "
                              "1.198793496548523e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 8       2 5    "
                              "5.777886299654152e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 8       2 6    "
                              "1.246636358465198e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 8       2 7   "
                              "-8.182168002205954e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 8       2 8   "
                              "-8.152060972764459e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 8       2 9    "
                              "6.397272639737098e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 8      2 10   "
                              "-9.673090992867567e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8       2 2    "
                              "2.299637577010168e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8       2 3   "
                              "-8.386188606718816e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8       2 4   "
                              "-3.546585198760869e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8       2 5   "
                              "-1.594447283167996e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8       2 6   "
                              "-3.624451798408940e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8       2 7    "
                              "2.393694254451009e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8       2 8    "
                              "2.380772494023351e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8       2 9   "
                              "-1.873508452455056e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 8      2 10    "
                              "2.792924338397144e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 1   "
                              "-1.482398297344698e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 2    "
                              "4.334481872461597e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 3   "
                              "-1.576891164176070e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 4   "
                              "-6.683123660819741e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 5   "
                              "-2.955970298514095e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 6   "
                              "-6.802638916863286e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 7    "
                              "4.498850162058932e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 8    "
                              "4.472979702640088e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8       2 9   "
                              "-3.521929716656321e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 8      2 10    "
                              "5.233411582492333e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8       2 2    "
                              "2.742839545884807e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8       2 3   "
                              "-1.007254829333757e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8       2 4   "
                              "-4.238545057657451e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8       2 5   "
                              "-1.959727484419271e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8       2 6   "
                              "-4.359343171640758e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8       2 7    "
                              "2.870595717144795e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8       2 8    "
                              "2.858295831875154e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8       2 9   "
                              "-2.245129471809810e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 8      2 10    "
                              "3.366992515594729e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 8       2 3    "
                              "2.785479656651575e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 8       2 4    "
                              "1.137894864898774e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 8       2 5    "
                              "6.237387858775510e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 8       2 6    "
                              "1.222199276527281e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 8       2 7   "
                              "-7.908438434585076e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 8       2 8   "
                              "-7.921312454946876e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 8       2 9    "
                              "6.161858483127646e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 8      2 10   "
                              "-9.591808561875573e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 8       2 6    "
                              "5.074230704757889e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 8       2 7   "
                              "-3.355787621948815e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 8       2 8   "
                              "-3.319416623364573e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 8       2 9    "
                              "2.636607645253481e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 8      2 10   "
                              "-3.974403634290948e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 9       2 3   "
                              "-1.482495648159440e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 9       2 4    "
                              "4.517805884059677e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 9       2 5   "
                              "-1.366111317915933e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 9       2 6    "
                              "1.520722898415504e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 9       2 7   "
                              "-6.061114169503505e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 9       2 8    "
                              "6.397272639737097e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 9       2 9    "
                              "6.581877662889447e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 0       2 9      2 10    "
                              "2.769873474544294e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9       2 2    "
                              "1.497673462499832e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9       2 3    "
                              "4.376268190834142e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9       2 4   "
                              "-1.333820104755788e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9       2 5    "
                              "4.025856583911213e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9       2 6   "
                              "-4.493517296103258e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9       2 7    "
                              "1.772110467094541e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9       2 8   "
                              "-1.873508452455056e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9       2 9   "
                              "-1.933943941591720e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 2       2 9      2 10   "
                              "-8.068081874153684e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9       2 2    "
                              "2.821190975353544e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9       2 3    "
                              "8.240513849039941e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9       2 4   "
                              "-2.512010284802130e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9       2 5    "
                              "7.577926689386553e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9       2 6   "
                              "-8.464039591527386e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9       2 7    "
                              "3.329629255617322e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9       2 8   "
                              "-3.521929716656321e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9       2 9   "
                              "-3.638031331079234e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4       2 9      2 10   "
                              "-1.514802333492648e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9       2 2    "
                              "1.789233425378196e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9       2 3    "
                              "5.219917234123755e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9       2 4   "
                              "-1.593432574526108e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9       2 5    "
                              "4.806215534947036e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9       2 6   "
                              "-5.364012338293060e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9       2 7    "
                              "2.122274753797954e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9       2 8   "
                              "-2.245129471809810e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9       2 9   "
                              "-2.312903290591723e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6       2 9      2 10   "
                              "-9.689036069479789e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 9       2 3   "
                              "-1.395412567391568e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 9       2 4    "
                              "4.282577880577520e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 9       2 5   "
                              "-1.291638838655702e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 9       2 6    "
                              "1.436182314881395e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 9       2 7   "
                              "-5.825026979372000e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 9       2 8    "
                              "6.161858483127646e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 9       2 9    "
                              "6.278283290908751e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 8       2 9      2 10    "
                              "2.696606804509703e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 9       2 7   "
                              "-2.542766261503350e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 9       2 8    "
                              "2.636607645253481e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 9       2 9    "
                              "2.743726361158304e+01 \n";
    integralFileTwoBodyFAD << "       1 4      1 10       2 9      2 10    "
                              "1.140490862655844e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0      2 10       2 4   "
                              "-5.492070086242435e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0      2 10       2 5    "
                              "2.130802009503947e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 0      2 10       2 6   "
                              "-2.250794069603903e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0      2 10       2 7    "
                              "1.094561935310272e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0      2 10       2 8   "
                              "-9.673090992867567e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0      2 10       2 9    "
                              "2.769873474544295e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 0      2 10      2 10   "
                              "-1.080485149673152e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10       2 2   "
                              "-1.525628551227822e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10       2 3    "
                              "2.393366749158278e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10       2 4    "
                              "1.609440180592267e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10       2 5   "
                              "-6.350926836375680e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10       2 6    "
                              "6.571034680643295e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10       2 7   "
                              "-3.216256790135352e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10       2 8    "
                              "2.792924338397145e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10       2 9   "
                              "-8.068081874153684e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 2      2 10      2 10    "
                              "3.145579135103967e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10       2 2   "
                              "-2.871439969097597e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10       2 3    "
                              "4.488331227771438e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10       2 4    "
                              "3.026180852572170e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10       2 5   "
                              "-1.198518721277427e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10       2 6    "
                              "1.234485106730705e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10       2 7   "
                              "-6.050833727024074e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10       2 8    "
                              "5.233411582492337e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10       2 9   "
                              "-1.514802333492648e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 4      2 10      2 10    "
                              "5.905590478277218e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10       2 2   "
                              "-1.825769883946854e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10       2 3    "
                              "2.886987091531497e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10       2 4    "
                              "1.929757887422544e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10       2 5   "
                              "-7.551003996121392e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10       2 6    "
                              "7.891806030829301e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10       2 7   "
                              "-3.850153106984002e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10       2 8    "
                              "3.366992515594731e+02 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10       2 9   "
                              "-9.689036069479791e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 6      2 10      2 10    "
                              "3.780849035335164e+03 \n";
    integralFileTwoBodyFAD << "       1 4       1 8      2 10       2 4   "
                              "-5.302150443587372e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8      2 10       2 5    "
                              "1.971049907245200e+00 \n";
    integralFileTwoBodyFAD << "       1 4       1 8      2 10       2 6   "
                              "-2.190682304157160e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8      2 10       2 7    "
                              "1.048492818872038e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8      2 10       2 8   "
                              "-9.591808561875575e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8      2 10       2 9    "
                              "2.696606804509704e+01 \n";
    integralFileTwoBodyFAD << "       1 4       1 8      2 10      2 10   "
                              "-1.055920326425384e+03 \n";
    integralFileTwoBodyFAD << "       1 4      1 10      2 10       2 8   "
                              "-3.974403634290948e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10      2 10       2 9    "
                              "1.140490862655845e+00 \n";
    integralFileTwoBodyFAD << "       1 4      1 10      2 10      2 10   "
                              "-4.409889931930822e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 0       2 1    "
                              "8.823396009909219e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 0       2 2   "
                              "-6.743992800579893e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 0       2 0   "
                              "-3.316175259432389e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 0       2 1    "
                              "4.594621739182922e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 0       2 2   "
                              "-3.417105912761267e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 0       2 3   "
                              "-2.350738657084364e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 0       2 0   "
                              "-5.354025400196306e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 0       2 1    "
                              "7.198080576474733e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 0       2 2   "
                              "-5.331586085740724e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 0       2 3   "
                              "-3.551134474409608e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 0       2 0    "
                              "3.586983501749820e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 0       2 1   "
                              "-5.596348316474126e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 0       2 2    "
                              "4.223994946478371e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 0       2 3    "
                              "3.089803528482331e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 0       2 1    "
                              "1.110179599027563e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 1       2 0    "
                              "8.823396009909217e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 1       2 2   "
                              "-1.249099522197051e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 1       2 3   "
                              "-1.148745155732086e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 1       2 4   "
                              "-1.822899804609036e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 1       2 0    "
                              "4.594621739182922e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 1       2 2   "
                              "-6.498731144299226e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 1       2 3   "
                              "-5.820462306712188e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 1       2 4   "
                              "-8.574476810455801e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 1       2 5    "
                              "1.193383332023027e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 1       2 6   "
                              "-1.978363485299689e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 1       2 0    "
                              "7.198080576474733e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 1       2 2   "
                              "-1.017926032868012e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 1       2 3   "
                              "-9.080690101768228e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 1       2 4   "
                              "-1.315652114516399e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 1       2 5    "
                              "1.843850815936944e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 1       2 6   "
                              "-3.119554296238144e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 1       2 7    "
                              "1.089193636956599e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 1       2 8   "
                              "-1.384831461383535e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 1       2 0   "
                              "-5.596348316474127e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 1       2 1   "
                              "-1.836616689658290e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 1       2 2    "
                              "7.918326825421234e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 1       2 3    "
                              "7.193974856486150e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 1       2 4    "
                              "1.092206969244397e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 1       2 5   "
                              "-1.468160434437427e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 1       2 6    "
                              "2.389380046566101e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 1       2 8    "
                              "1.081104317671499e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 1       2 0    "
                              "1.110179599027563e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 1       2 2   "
                              "-1.572565038600372e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 1       2 3   "
                              "-1.588406858635852e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 2       2 0   "
                              "-6.743992800579893e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 2       2 1   "
                              "-1.249099522197051e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 2       2 2    "
                              "1.588634273425815e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 2       2 3    "
                              "1.509761598759491e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 2       2 4    "
                              "1.845337029158175e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 0   "
                              "-3.417105912761269e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 1   "
                              "-6.498731144299226e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 2    "
                              "4.978911547778019e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 3    "
                              "7.847674701878742e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 4    "
                              "9.363634812346245e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 5    "
                              "3.200745813643807e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 7    "
                              "3.358013897913122e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 8    "
                              "2.576579352025349e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2       2 9    "
                              "1.678859903465051e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 2      2 10   "
                              "-1.711175652333371e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 0   "
                              "-5.331586085740724e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 1   "
                              "-1.017926032868012e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 2    "
                              "7.037620859735704e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 3    "
                              "1.228980012923533e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 4    "
                              "1.460999000485207e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 5    "
                              "5.339230508476397e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 7    "
                              "5.270721006979482e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 8    "
                              "4.046669106670033e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2       2 9    "
                              "2.634358780813928e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 2      2 10   "
                              "-2.681576095260986e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 0    "
                              "4.223994946478371e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 1    "
                              "7.918326825421234e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 2   "
                              "-8.209974031941123e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 3   "
                              "-9.565319590112397e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 4   "
                              "-1.156118632219219e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 5   "
                              "-3.349919430076855e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 7   "
                              "-4.080648406830757e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 8   "
                              "-3.119557424450956e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2       2 9   "
                              "-2.038515382809466e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 2      2 10    "
                              "2.084417658589889e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 2       2 1   "
                              "-1.572565038600372e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 2       2 3    "
                              "1.901626899321122e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 2       2 4    "
                              "2.525909890768939e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 3       2 1   "
                              "-1.148745155732085e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 3       2 2    "
                              "1.509761598759491e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 3       2 3   "
                              "-3.000836027715836e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 3       2 4    "
                              "1.829035662226727e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 3       2 5   "
                              "-7.801003823052968e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 3       2 6   "
                              "-1.362561801522244e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 3       2 7    "
                              "2.760620515744144e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 3       2 8   "
                              "-1.825960811296089e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 0   "
                              "-2.350738657084360e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 1   "
                              "-5.820462306712191e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 2    "
                              "7.847674701878742e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 3   "
                              "-2.034217450612956e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 4    "
                              "9.498941090852404e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 5   "
                              "-3.863105522499197e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 6   "
                              "-6.827676857536898e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 7    "
                              "1.452537613796562e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 8   "
                              "-9.415227099173276e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3       2 9    "
                              "4.904632855030777e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 3      2 10    "
                              "2.691417579637120e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 0   "
                              "-3.551134474409590e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 1   "
                              "-9.080690101768234e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 2    "
                              "1.228980012923533e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 3   "
                              "-3.296807542592913e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 4    "
                              "1.487333306297969e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 5   "
                              "-6.008135788476922e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 6   "
                              "-1.061728191681495e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 7    "
                              "2.277864278398907e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 8   "
                              "-1.473526625578310e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3       2 9    "
                              "7.690058997159826e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 3      2 10    "
                              "4.194139466786696e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 0    "
                              "3.089803528482317e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 1    "
                              "7.193974856486152e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 2   "
                              "-9.565319590112397e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 3    "
                              "2.171104921565680e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 4   "
                              "-1.158248083248140e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 5    "
                              "4.841163461759432e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 6    "
                              "8.459141036347970e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 7   "
                              "-1.757941117076523e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 8    "
                              "1.153836792718763e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3       2 9   "
                              "-5.943802001611473e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 3      2 10   "
                              "-3.325931811511963e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 3       2 1   "
                              "-1.588406858635852e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 3       2 2    "
                              "1.901626899321123e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 3       2 4    "
                              "2.307273346932320e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 3       2 5   "
                              "-1.183947176249873e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 3       2 6   "
                              "-1.820141439304294e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 1   "
                              "-1.822899804609030e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 2    "
                              "1.845337029158175e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 3    "
                              "1.829035662226727e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 4    "
                              "3.226385831330303e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 5   "
                              "-1.585052114112287e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 6   "
                              "-6.286568347067822e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 7   "
                              "-7.436477558802694e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 8   "
                              "-7.564426310747169e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4       2 9   "
                              "-2.852205906660334e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 4      2 10    "
                              "3.495015082827114e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 1   "
                              "-8.574476810455772e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 2    "
                              "9.363634812346251e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 3    "
                              "9.498941090852404e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 4    "
                              "1.625842713109470e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 5   "
                              "-8.231351574909822e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 6   "
                              "-3.219045746185515e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 7   "
                              "-3.964321159021027e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 8   "
                              "-3.975455982742326e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4       2 9   "
                              "-1.495263735738211e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 4      2 10    "
                              "1.806368196341459e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 1   "
                              "-1.315652114516397e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 2    "
                              "1.460999000485207e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 3    "
                              "1.487333306297969e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 4    "
                              "2.534058924154749e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 5   "
                              "-1.288728509165082e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 6   "
                              "-5.028320235850737e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 7   "
                              "-6.238133268191687e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 8   "
                              "-6.242237908992061e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4       2 9   "
                              "-2.345409585059690e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 4      2 10    "
                              "2.826352610167659e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 1    "
                              "1.092206969244393e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 2   "
                              "-1.156118632219221e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 3   "
                              "-1.158248083248140e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 4   "
                              "-2.014944216893012e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 5    "
                              "1.003491365462847e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 6    "
                              "3.952745899522264e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 7    "
                              "4.781871010178232e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 8    "
                              "4.828240276279016e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4       2 9    "
                              "1.815838556879358e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 4      2 10   "
                              "-2.208292144964529e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 4       2 2    "
                              "2.525909890768939e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 4       2 3    "
                              "2.307273346932320e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 4       2 4    "
                              "4.516657380965690e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 4       2 5   "
                              "-1.990417643706983e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 4       2 6   "
                              "-8.241155019424012e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 5       2 3   "
                              "-7.801003823052961e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 5       2 4   "
                              "-1.585052114112287e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 5       2 5   "
                              "-1.434643881367362e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 5       2 6    "
                              "1.940177984279219e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 5       2 7    "
                              "3.685733728801261e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 5       2 8   "
                              "-3.881736394290465e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 5       2 9    "
                              "8.626285889977325e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 5      2 10   "
                              "-1.331368796000002e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 1    "
                              "1.193383332022999e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 2    "
                              "3.200745813643785e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 3   "
                              "-3.863105522499191e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 4   "
                              "-8.231351574909822e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 5   "
                              "-7.520633571223505e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 6    "
                              "1.004750511298486e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 7    "
                              "1.918985185655001e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 8   "
                              "-1.805567291700113e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5       2 9    "
                              "4.513174695465784e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 5      2 10   "
                              "-7.108706853343747e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 1    "
                              "1.843850815936916e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 2    "
                              "5.339230508476397e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 3   "
                              "-6.008135788476922e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 4   "
                              "-1.288728509165082e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 5   "
                              "-1.179195059392546e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 6    "
                              "1.572402092689695e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 7    "
                              "3.004214286091338e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 8   "
                              "-2.762931093050269e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5       2 9    "
                              "7.072332145784728e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 5      2 10   "
                              "-1.118572195141371e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 1   "
                              "-1.468160434437413e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 2   "
                              "-3.349919430076834e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 3    "
                              "4.841163461759426e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 4    "
                              "1.003491365462848e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 5    "
                              "9.124093220817824e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 6   "
                              "-1.226713304451232e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 7   "
                              "-2.333108475191312e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 8    "
                              "2.311902370230600e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5       2 9   "
                              "-5.478109977465269e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 5      2 10    "
                              "8.558681798928372e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 5       2 3   "
                              "-1.183947176249875e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 5       2 4   "
                              "-1.990417643706982e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 5       2 5   "
                              "-1.745176364199177e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 5       2 6    "
                              "2.460688425725433e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 5       2 7    "
                              "4.460267096941003e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 5       2 9    "
                              "1.043158182951087e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 6       2 3   "
                              "-1.362561801522245e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 6       2 4   "
                              "-6.286568347067822e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 6       2 5    "
                              "1.940177984279218e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 6       2 6    "
                              "2.631392050094647e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 6       2 7   "
                              "-8.382115338489656e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 6       2 8   "
                              "-7.991987013406887e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 6       2 9   "
                              "-9.587835720001364e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 6      2 10    "
                              "1.437706639836916e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6       2 1   "
                              "-1.978363485299689e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6       2 3   "
                              "-6.827676857536898e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6       2 4   "
                              "-3.219045746185515e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6       2 5    "
                              "1.004750511298486e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6       2 6    "
                              "1.356266944253072e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6       2 7   "
                              "-4.348436915246994e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6       2 8   "
                              "-4.072545088132770e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6       2 9   "
                              "-5.036452631835036e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 6      2 10    "
                              "7.379221537115295e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6       2 1   "
                              "-3.119554296238103e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6       2 3   "
                              "-1.061728191681495e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6       2 4   "
                              "-5.028320235850737e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6       2 5    "
                              "1.572402092689695e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6       2 6    "
                              "2.120971270444157e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6       2 7   "
                              "-6.806227284655689e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6       2 8   "
                              "-6.353848730438436e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6       2 9   "
                              "-7.901502432058248e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 6      2 10    "
                              "1.153066208706179e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6       2 1    "
                              "2.389380046566102e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6       2 3    "
                              "8.459141036347962e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6       2 4    "
                              "3.952745899522264e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6       2 5   "
                              "-1.226713304451232e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6       2 6   "
                              "-1.660072126397794e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6       2 7    "
                              "5.302471945312025e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6       2 8    "
                              "5.008491166851751e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6       2 9    "
                              "6.108645948415147e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 6      2 10   "
                              "-9.048897347514762e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 6       2 3   "
                              "-1.820141439304296e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 6       2 4   "
                              "-8.241155019424008e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 6       2 5    "
                              "2.460688425725433e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 6       2 6    "
                              "3.393520885342537e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 6       2 7   "
                              "-1.049786435114059e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 6       2 8   "
                              "-1.041815763577076e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 6       2 9   "
                              "-1.177377020620830e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 6      2 10    "
                              "1.866941760448757e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 7       2 3    "
                              "2.760620515744144e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 7       2 4   "
                              "-7.436477558802687e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 7       2 5    "
                              "3.685733728801262e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 7       2 6   "
                              "-8.382115338489656e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 7       2 7   "
                              "-3.296506992663207e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 7       2 8    "
                              "5.212488665443412e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 7       2 9    "
                              "3.858025883183871e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 7      2 10   "
                              "-6.943714719025532e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7       2 2    "
                              "3.358013897913122e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7       2 3    "
                              "1.452537613796562e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7       2 4   "
                              "-3.964321159021035e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7       2 5    "
                              "1.918985185655001e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7       2 6   "
                              "-4.348436915246994e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7       2 7   "
                              "-1.718113719975473e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7       2 8    "
                              "2.686979812925845e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7       2 9    "
                              "1.988655440864239e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 7      2 10   "
                              "-3.608016257964326e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 1    "
                              "1.089193636956610e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 2    "
                              "5.270721006979511e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 3    "
                              "2.277864278398907e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 4   "
                              "-6.238133268191687e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 5    "
                              "3.004214286091338e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 6   "
                              "-6.806227284655689e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 7   "
                              "-2.691259231743191e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 8    "
                              "4.201183475929211e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7       2 9    "
                              "3.107888834075448e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 7      2 10   "
                              "-5.650171852750259e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7       2 2   "
                              "-4.080648406830750e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7       2 3   "
                              "-1.757941117076524e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7       2 4    "
                              "4.781871010178232e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7       2 5   "
                              "-2.333108475191312e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7       2 6    "
                              "5.302471945312025e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7       2 7    "
                              "2.090162750628617e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7       2 8   "
                              "-3.286902251655388e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7       2 9   "
                              "-2.428842699671646e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 7      2 10    "
                              "4.398595092979114e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 7       2 5    "
                              "4.460267096940998e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 7       2 6   "
                              "-1.049786435114059e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 7       2 7   "
                              "-4.070465056235971e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 7       2 8    "
                              "6.645736996135043e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 7       2 9    "
                              "4.774751388295480e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 8       2 3   "
                              "-1.825960811296086e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 8       2 4   "
                              "-7.564426310747154e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 8       2 5   "
                              "-3.881736394290473e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 8       2 6   "
                              "-7.991987013406886e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 8       2 7    "
                              "5.212488665443412e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 8       2 8    "
                              "5.204120180804815e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 8       2 9   "
                              "-4.070085339804572e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 8      2 10    "
                              "6.237789957087999e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8       2 2    "
                              "2.576579352025348e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8       2 3   "
                              "-9.415227099173276e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8       2 4   "
                              "-3.975455982742323e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8       2 5   "
                              "-1.805567291700107e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8       2 6   "
                              "-4.072545088132770e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8       2 7    "
                              "2.686979812925845e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8       2 8    "
                              "2.673363584959902e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8       2 9   "
                              "-2.102601543823059e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 8      2 10    "
                              "3.141131878579666e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 1   "
                              "-1.384831461383556e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 2    "
                              "4.046669106670033e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 3   "
                              "-1.473526625578302e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 4   "
                              "-6.242237908992047e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 5   "
                              "-2.762931093050274e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 6   "
                              "-6.353848730438438e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 7    "
                              "4.201183475929212e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 8    "
                              "4.177659848343106e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8       2 9   "
                              "-3.288556693407418e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 8      2 10    "
                              "4.887741685646961e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 1    "
                              "1.081104317671478e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 2   "
                              "-3.119557424450954e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 3    "
                              "1.153836792718762e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 4    "
                              "4.828240276279022e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 5    "
                              "2.311902370230600e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 6    "
                              "5.008491166851751e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 7   "
                              "-3.286902251655388e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 8   "
                              "-3.276419930057813e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8       2 9    "
                              "2.568935552594200e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 8      2 10   "
                              "-3.881011585810597e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 8       2 6   "
                              "-1.041815763577076e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 8       2 7    "
                              "6.645736996135042e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 8       2 8    "
                              "6.718945757846798e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 8       2 9   "
                              "-5.144945807610153e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 8      2 10    "
                              "8.157453662324443e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 9       2 4   "
                              "-2.852205906660329e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 9       2 5    "
                              "8.626285889977346e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 9       2 6   "
                              "-9.587835720001371e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 9       2 7    "
                              "3.858025883183871e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 9       2 8   "
                              "-4.070085339804572e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 9       2 9   "
                              "-4.171402579764400e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 1       2 9      2 10   "
                              "-1.771431285625149e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9       2 2    "
                              "1.678859903465039e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9       2 3    "
                              "4.904632855030806e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9       2 4   "
                              "-1.495263735738216e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9       2 5    "
                              "4.513174695465778e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9       2 6   "
                              "-5.036452631835036e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9       2 7    "
                              "1.988655440864239e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9       2 8   "
                              "-2.102601543823059e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9       2 9   "
                              "-2.169100351846718e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 3       2 9      2 10   "
                              "-9.061377781098231e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9       2 2    "
                              "2.634358780813950e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9       2 3    "
                              "7.690058997159809e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9       2 4   "
                              "-2.345409585059683e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9       2 5    "
                              "7.072332145784728e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9       2 6   "
                              "-7.901502432058248e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9       2 7    "
                              "3.107888834075448e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9       2 8   "
                              "-3.288556693407418e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9       2 9   "
                              "-3.396122610559982e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5       2 9      2 10   "
                              "-1.414512481127354e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9       2 2   "
                              "-2.038515382809473e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9       2 3   "
                              "-5.943802001611459e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9       2 4    "
                              "1.815838556879358e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9       2 5   "
                              "-5.478109977465269e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9       2 6    "
                              "6.108645948415145e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9       2 7   "
                              "-2.428842699671645e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9       2 8    "
                              "2.568935552594200e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9       2 9    "
                              "2.641087236144967e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7       2 9      2 10    "
                              "1.111693254863085e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 9       2 5    "
                              "1.043158182951086e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 9       2 6   "
                              "-1.177377020620829e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 9       2 7    "
                              "4.774751388295484e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 9       2 8   "
                              "-5.144945807610153e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 9       2 9   "
                              "-5.157927070990900e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 9       2 9      2 10   "
                              "-2.266214214577093e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1      2 10       2 4    "
                              "3.495015082827114e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1      2 10       2 5   "
                              "-1.331368796000000e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1      2 10       2 6    "
                              "1.437706639836915e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1      2 10       2 7   "
                              "-6.943714719025532e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 1      2 10       2 8    "
                              "6.237789957087999e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1      2 10       2 9   "
                              "-1.771431285625149e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 1      2 10      2 10    "
                              "6.917099271935308e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10       2 2   "
                              "-1.711175652333362e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10       2 3    "
                              "2.691417579637126e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10       2 4    "
                              "1.806368196341459e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10       2 5   "
                              "-7.108706853343762e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10       2 6    "
                              "7.379221537115293e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10       2 7   "
                              "-3.608016257964327e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10       2 8    "
                              "3.141131878579664e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10       2 9   "
                              "-9.061377781098231e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 3      2 10      2 10    "
                              "3.533710787834432e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10       2 2   "
                              "-2.681576095260979e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10       2 3    "
                              "4.194139466786687e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10       2 4    "
                              "2.826352610167657e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10       2 5   "
                              "-1.118572195141369e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10       2 6    "
                              "1.153066208706179e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10       2 7   "
                              "-5.650171852750260e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10       2 8    "
                              "4.887741685646959e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10       2 9   "
                              "-1.414512481127354e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 5      2 10      2 10    "
                              "5.515661074536773e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10       2 2    "
                              "2.084417658589892e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10       2 3   "
                              "-3.325931811511963e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10       2 4   "
                              "-2.208292144964529e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10       2 5    "
                              "8.558681798928356e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10       2 6   "
                              "-9.048897347514762e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10       2 7    "
                              "4.398595092979114e+01 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10       2 8   "
                              "-3.881011585810597e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10       2 9    "
                              "1.111693254863085e+02 \n";
    integralFileTwoBodyFAD << "       1 5       1 7      2 10      2 10   "
                              "-4.340608353924198e+03 \n";
    integralFileTwoBodyFAD << "       1 5       1 9      2 10       2 6    "
                              "1.866941760448757e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9      2 10       2 8    "
                              "8.157453662324443e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9      2 10       2 9   "
                              "-2.266214214577093e+00 \n";
    integralFileTwoBodyFAD << "       1 5       1 9      2 10      2 10    "
                              "8.964759125910373e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 0       2 1   "
                              "-2.726237669356156e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 0       2 2    "
                              "2.141076998089606e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 0       2 1    "
                              "9.847171376405946e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 0       2 2   "
                              "-7.571597286196535e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 0       2 0   "
                              "-3.372303250186400e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 0       2 1    "
                              "4.901474884730100e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 0       2 2   "
                              "-3.667946726863760e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 0       2 3   "
                              "-2.575222147295233e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 0       2 0   "
                              "-7.152302226954109e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 0       2 1    "
                              "1.031485594529488e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 0       2 2   "
                              "-7.711347896542597e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 0       2 3   "
                              "-5.386959103789573e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 0       2 0    "
                              "6.111356830875017e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 0       2 1   "
                              "-9.055431617871011e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 0       2 2    "
                              "6.792768343223963e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 0       2 3    "
                              "4.834720798096544e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 0       2 1    "
                              "1.172883015425946e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 0       2 2   "
                              "-9.142470620259777e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 1       2 0   "
                              "-2.726237669356156e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 1       2 2    "
                              "3.860781707625213e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 1       2 3    "
                              "3.644759010542844e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 1       2 0    "
                              "9.847171376405947e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 1       2 2   "
                              "-1.394109019314518e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 1       2 3   "
                              "-1.289598697485621e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 1       2 4   "
                              "-2.054573374740645e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 1       2 0    "
                              "4.901474884730099e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 1       2 2   "
                              "-6.933484505246231e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 1       2 3   "
                              "-6.247104701786767e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 1       2 4   "
                              "-9.289958677014988e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 1       2 5    "
                              "1.272510768566129e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 1       2 6   "
                              "-2.106840597584706e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 0    "
                              "1.031485594529488e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 1    "
                              "1.792725440805916e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 2   "
                              "-1.459068171657053e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 3   "
                              "-1.313361667380778e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 4   "
                              "-1.948172777994673e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 5    "
                              "2.673727298675128e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 6   "
                              "-4.437617660052866e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 7    "
                              "1.551477429913042e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 1       2 8   "
                              "-1.988457274011611e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 0   "
                              "-9.055431617871011e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 1   "
                              "-2.067448858236699e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 2    "
                              "1.281057154462392e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 3    "
                              "1.156926370480235e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 4    "
                              "1.732519626425011e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 5   "
                              "-2.361289543442362e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 6    "
                              "3.882856635079146e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 7   "
                              "-1.359032754209603e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 1       2 8    "
                              "1.746843528431626e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 1       2 0    "
                              "1.172883015425946e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 1       2 1    "
                              "1.008383095100735e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 1       2 2   "
                              "-1.660150637030747e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 1       2 3   "
                              "-1.555830872314816e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 1       2 4   "
                              "-2.434022040650357e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 2       2 0    "
                              "2.141076998089606e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 2       2 1    "
                              "3.860781707625213e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 2       2 3   "
                              "-4.668008245825254e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 2       2 4   "
                              "-5.840263653569882e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 2       2 0   "
                              "-7.571597286196535e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 2       2 1   "
                              "-1.394109019314518e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 2       2 2    "
                              "1.929094579319154e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 2       2 3    "
                              "1.685115800648999e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 2       2 4    "
                              "2.070443939231001e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 0   "
                              "-3.667946726863763e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 1   "
                              "-6.933484505246231e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 2    "
                              "6.093738515021593e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 3    "
                              "8.373557364342490e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 4    "
                              "1.004477670517755e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 5    "
                              "3.253023759749116e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 7    "
                              "3.580764339138637e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 8    "
                              "2.742839545884832e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2       2 9    "
                              "1.789233425378168e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 2      2 10   "
                              "-1.825769883946858e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 0   "
                              "-7.711347896542597e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 1   "
                              "-1.459068171657053e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 2    "
                              "1.255830495875296e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 3    "
                              "1.762062211703427e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 4    "
                              "2.111884212475056e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 5    "
                              "6.924399637986721e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 7    "
                              "7.537250142642635e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 8    "
                              "5.774732962362257e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2       2 9    "
                              "3.766287808970467e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 2      2 10   "
                              "-3.842278465414013e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 0    "
                              "6.792768343223963e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 1    "
                              "1.281057154462392e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 2   "
                              "-1.182443949774317e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 3   "
                              "-1.547259245000407e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 4   "
                              "-1.860006153611544e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 5   "
                              "-5.821056516133519e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 7   "
                              "-6.610698557937501e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 8   "
                              "-5.061222867872859e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2       2 9   "
                              "-3.303255145773667e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 2      2 10    "
                              "3.372852130492381e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 2       2 0   "
                              "-9.142470620259777e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 2       2 1   "
                              "-1.660150637030747e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 2       2 2    "
                              "2.719105181960331e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 2       2 3    "
                              "2.006180858615087e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 2       2 4    "
                              "2.492695068473067e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 3       2 1    "
                              "3.644759010542840e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 3       2 2   "
                              "-4.668008245825254e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 3       2 4   "
                              "-5.658183131020964e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 3       2 5    "
                              "2.541196880784557e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 3       2 6    "
                              "4.303859263961107e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 1   "
                              "-1.289598697485620e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 2    "
                              "1.685115800648999e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 3   "
                              "-3.125387959776451e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 4    "
                              "2.041680613739775e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 5   "
                              "-8.815444843999501e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 6   "
                              "-1.527330805382476e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 7    "
                              "3.071950060034313e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 8   "
                              "-2.043744512464690e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 3       2 9    "
                              "1.040082899463002e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 0   "
                              "-2.575222147295244e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 1   "
                              "-6.247104701786767e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 2    "
                              "8.373557364342490e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 3   "
                              "-2.058523508936851e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 4    "
                              "1.013686331149331e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 5   "
                              "-4.172003642850805e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 6   "
                              "-7.327405851298587e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 7    "
                              "1.544951744941340e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 8   "
                              "-1.007254829333753e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3       2 9    "
                              "5.219917234123790e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 3      2 10    "
                              "2.886987091531478e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 0   "
                              "-5.386959103789601e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 1   "
                              "-1.313361667380778e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 2    "
                              "1.762062211703427e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 3   "
                              "-4.370206185986421e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 4    "
                              "2.133055333698510e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 5   "
                              "-8.762977869065867e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 6   "
                              "-1.539980399813393e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 7    "
                              "3.252395435257811e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 8   "
                              "-2.118969104465794e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3       2 9    "
                              "1.098827415807879e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 3      2 10    "
                              "6.069100412154372e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 0    "
                              "4.834720798096551e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 1    "
                              "1.156926370480236e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 2   "
                              "-1.547259245000407e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 3    "
                              "3.721939784918222e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 4   "
                              "-1.873227694276366e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 5    "
                              "7.741903802694505e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 6    "
                              "1.358448459682019e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 7   "
                              "-2.851818295367525e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 8    "
                              "1.862452156449607e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3       2 9   "
                              "-9.637227239216637e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 3      2 10   "
                              "-5.348072594774349e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 3       2 4    "
                              "1.167903155640119e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 1   "
                              "-1.555830872314816e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 2    "
                              "2.006180858615087e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 3   "
                              "-3.131303957664550e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 4    "
                              "2.430724079773960e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 5   "
                              "-1.080761725485800e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 6   "
                              "-1.818513046857066e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 7    "
                              "3.618890254776641e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 8   "
                              "-2.459244785532948e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 3       2 9    "
                              "1.228450458146215e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 4       2 2   "
                              "-5.840263653569882e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 4       2 3   "
                              "-5.658183131020966e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 4       2 4   "
                              "-1.027814612081334e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 4       2 5    "
                              "4.899019466432279e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 4       2 6    "
                              "1.967239206932404e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 4       2 7    "
                              "2.269065133179011e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 4       2 8    "
                              "2.336647615369954e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 4      2 10   "
                              "-1.086434005727397e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 1   "
                              "-2.054573374740664e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 2    "
                              "2.070443939231001e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 3    "
                              "2.041680613739775e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 4    "
                              "3.625385721839491e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 5   "
                              "-1.768902579284091e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 6   "
                              "-7.035113844255113e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 7   "
                              "-8.283479859434099e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 8   "
                              "-8.442497964066138e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4       2 9   "
                              "-3.180643948960411e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 4      2 10    "
                              "3.905523455266483e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 1   "
                              "-9.289958677014960e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 2    "
                              "1.004477670517754e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 3    "
                              "1.013686331149331e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 4    "
                              "1.746846391965047e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 5   "
                              "-8.782855812483540e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 6   "
                              "-3.444658859938158e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 7   "
                              "-4.215314041503886e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 8   "
                              "-4.238545057657445e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4       2 9   "
                              "-1.593432574526110e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 4      2 10    "
                              "1.929757887422542e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 1   "
                              "-1.948172777994688e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 2    "
                              "2.111884212475056e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 3    "
                              "2.133055333698510e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 4    "
                              "3.671750516312780e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 5   "
                              "-1.848145744397441e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 6   "
                              "-7.244790812122642e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 7   "
                              "-8.877390224162222e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 8   "
                              "-8.922444358257381e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4       2 9   "
                              "-3.354025486260510e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 4      2 10    "
                              "4.059993790292226e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 1    "
                              "1.732519626425011e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 2   "
                              "-1.860006153611544e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 3   "
                              "-1.873227694276366e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 4   "
                              "-3.236647582216826e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 5    "
                              "1.623031378217881e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 6    "
                              "6.373507062916848e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 7    "
                              "7.771610193093105e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 8    "
                              "7.824021992367236e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4       2 9    "
                              "2.942090964056192e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 4      2 10   "
                              "-3.567575078941882e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 4       2 3    "
                              "1.167903155640119e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 4       2 5   "
                              "-1.010903912981590e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 1   "
                              "-2.434022040650364e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 2    "
                              "2.492695068473067e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 3    "
                              "2.430724079773960e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 4    "
                              "4.378540156496530e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 5   "
                              "-2.103620981509874e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 6   "
                              "-8.407319943652750e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 7   "
                              "-9.883299497889379e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 8   "
                              "-1.012354241756133e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4       2 9   "
                              "-3.787129376070258e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 4      2 10    "
                              "4.660740009072176e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 5       2 3    "
                              "2.541196880784560e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 5       2 4    "
                              "4.899019466432279e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 5       2 5    "
                              "4.395434511131091e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 5       2 6   "
                              "-6.012988730453160e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 5       2 7   "
                              "-1.130060948089802e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 5       2 8    "
                              "1.268994946619717e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 5       2 9   "
                              "-2.641232153392844e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 5       2 3   "
                              "-8.815444843999501e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 5       2 4   "
                              "-1.768902579284091e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 5       2 5   "
                              "-1.597955717538821e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 5       2 6    "
                              "2.166537252638214e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 5       2 7    "
                              "4.103989810285587e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 5       2 8   "
                              "-4.382514738828246e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 5       2 9    "
                              "9.607676452642700e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 5      2 10   "
                              "-1.480250734659055e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 1    "
                              "1.272510768566129e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 2    "
                              "3.253023759749116e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 3   "
                              "-4.172003642850810e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 4   "
                              "-8.782855812483539e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 5   "
                              "-8.008799078254622e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 6    "
                              "1.072723553863097e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 7    "
                              "2.044489550854506e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 8   "
                              "-1.959727484419270e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5       2 9    "
                              "4.806215534947030e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 5      2 10   "
                              "-7.551003996121406e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 1    "
                              "2.673727298675071e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 2    "
                              "6.924399637986650e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 3   "
                              "-8.762977869065890e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 4   "
                              "-1.848145744397441e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 5   "
                              "-1.685832711547716e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 6    "
                              "2.257065916510983e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 7    "
                              "4.302803370849783e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 8   "
                              "-4.107677408558210e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5       2 9    "
                              "1.011625124175042e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 5      2 10   "
                              "-1.590427668219140e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 1   "
                              "-2.361289543442359e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 2   "
                              "-5.821056516133477e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 3    "
                              "7.741903802694516e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 4    "
                              "1.623031378217881e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 5    "
                              "1.478768439004187e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 6   "
                              "-1.982828794678772e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 7   "
                              "-3.777338787717429e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 8    "
                              "3.658743194107578e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5       2 9   "
                              "-8.876029555580901e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 5      2 10    "
                              "1.391957884620224e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 5       2 4   "
                              "-1.010903912981590e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 5       2 6    "
                              "1.239387730666049e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 5       2 3   "
                              "-1.080761725485800e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 5       2 4   "
                              "-2.103620981509875e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 5       2 5   "
                              "-1.893254901229217e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 5       2 6    "
                              "2.579789406880491e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 5       2 7    "
                              "4.843737105543610e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 5       2 8   "
                              "-5.178464362429599e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 5       2 9    "
                              "1.135705037198194e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 5      2 10   "
                              "-1.757020148799413e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 6       2 3    "
                              "4.303859263961105e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 6       2 4    "
                              "1.967239206932402e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 6       2 5   "
                              "-6.012988730453159e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 6       2 6   "
                              "-8.192923728174320e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 6       2 7    "
                              "2.590026409798207e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 6       2 8    "
                              "2.501286650343349e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 6       2 9    "
                              "2.941297609416028e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 6      2 10   "
                              "-4.487282447436253e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 6       2 3   "
                              "-1.527330805382477e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 6       2 4   "
                              "-7.035113844255113e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 6       2 5    "
                              "2.166537252638214e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 6       2 6    "
                              "2.941473379694497e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 6       2 7   "
                              "-9.352934507908495e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 6       2 8   "
                              "-8.942954506089153e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 6       2 9   "
                              "-1.068556733649409e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 6      2 10    "
                              "1.607878057938480e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6       2 1   "
                              "-2.106840597584720e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6       2 3   "
                              "-7.327405851298590e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6       2 4   "
                              "-3.444658859938157e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6       2 5    "
                              "1.072723553863097e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6       2 6    "
                              "1.449533889084441e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6       2 7   "
                              "-4.639822205897414e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6       2 8   "
                              "-4.359343171640759e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6       2 9   "
                              "-5.364012338293058e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 6      2 10    "
                              "7.891806030829302e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6       2 1   "
                              "-4.437617660052811e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6       2 3   "
                              "-1.539980399813396e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6       2 4   "
                              "-7.244790812122642e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6       2 5    "
                              "2.257065916510983e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6       2 6    "
                              "3.049363897797529e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6       2 7   "
                              "-9.763169004791774e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6       2 8   "
                              "-9.167119239052861e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6       2 9   "
                              "-1.129153857177379e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 6      2 10    "
                              "1.659954508437191e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6       2 1    "
                              "3.882856635079158e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6       2 3    "
                              "1.358448459682021e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6       2 4    "
                              "6.373507062916848e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6       2 5   "
                              "-1.982828794678772e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6       2 6   "
                              "-2.680443857574112e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6       2 7    "
                              "8.574994994802636e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6       2 8    "
                              "8.069562405116351e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6       2 9    "
                              "9.902405437257963e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 6      2 10   "
                              "-1.459885354939898e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 6       2 5    "
                              "1.239387730666049e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 6       2 6    "
                              "1.685541059223176e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 6       2 3   "
                              "-1.818513046857068e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 6       2 4   "
                              "-8.407319943652750e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 6       2 5    "
                              "2.579789406880491e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 6       2 6    "
                              "3.510151522098442e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 6       2 7   "
                              "-1.111132929897306e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 6       2 8   "
                              "-1.065161737061374e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 6       2 9   "
                              "-1.269762453613440e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 6      2 10    "
                              "1.918822020574652e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 7       2 4    "
                              "2.269065133179007e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 7       2 5   "
                              "-1.130060948089803e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 7       2 6    "
                              "2.590026409798208e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 7       2 7    "
                              "1.014494595332452e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 7       2 8   "
                              "-1.619131998550263e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 7       2 9   "
                              "-1.191198100873518e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 7      2 10    "
                              "2.147959907634374e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 7       2 3    "
                              "3.071950060034313e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 7       2 4   "
                              "-8.283479859434113e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 7       2 5    "
                              "4.103989810285591e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 7       2 6   "
                              "-9.352934507908495e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 7       2 7   "
                              "-3.675033968435039e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 7       2 8    "
                              "5.823096020486471e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 7       2 9    "
                              "4.300596499617232e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 7      2 10   "
                              "-7.750490431838680e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7       2 2    "
                              "3.580764339138645e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7       2 3    "
                              "1.544951744941337e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7       2 4   "
                              "-4.215314041503891e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7       2 5    "
                              "2.044489550854507e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7       2 6   "
                              "-4.639822205897414e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7       2 7   "
                              "-1.831543874908211e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7       2 8    "
                              "2.870595717144795e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7       2 9    "
                              "2.122274753797954e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 7      2 10   "
                              "-3.850153106984001e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 1    "
                              "1.551477429913014e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 2    "
                              "7.537250142642479e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 3    "
                              "3.252395435257828e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 4   "
                              "-8.877390224162211e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 5    "
                              "4.302803370849783e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 6   "
                              "-9.763169004791775e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 7   "
                              "-3.854597791429401e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 8    "
                              "6.038952957782194e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7       2 9    "
                              "4.465101923264919e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 7      2 10   "
                              "-8.101896372023504e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 1   "
                              "-1.359032754209649e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 2   "
                              "-6.610698557937515e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 3   "
                              "-2.851818295367519e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 4    "
                              "7.771610193093099e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 5   "
                              "-3.777338787717432e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 6    "
                              "8.574994994802636e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 7    "
                              "3.383542788760048e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 8   "
                              "-5.308229091527654e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7       2 9   "
                              "-3.924516205341246e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 7      2 10    "
                              "7.114493677598679e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 7       2 7   "
                              "-2.096101407636399e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 7       2 3    "
                              "3.618890254776641e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 7       2 4   "
                              "-9.883299497889370e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 7       2 5    "
                              "4.843737105543613e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 7       2 6   "
                              "-1.111132929897306e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 7       2 7   "
                              "-4.359417594991877e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 7       2 8    "
                              "6.929917642364796e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 7       2 9    "
                              "5.083215609747403e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 7      2 10   "
                              "-9.231586122277200e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 8       2 4    "
                              "2.336647615369954e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 8       2 5    "
                              "1.268994946619717e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 8       2 6    "
                              "2.501286650343349e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 8       2 7   "
                              "-1.619131998550263e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 8       2 8   "
                              "-1.622133415357713e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 8       2 9    "
                              "1.261371567272619e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 8      2 10   "
                              "-1.959826459762755e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 8       2 3   "
                              "-2.043744512464701e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 8       2 4   "
                              "-8.442497964066153e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 8       2 5   "
                              "-4.382514738828245e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 8       2 6   "
                              "-8.942954506089153e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 8       2 7    "
                              "5.823096020486471e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 8       2 8    "
                              "5.818523303454713e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 8       2 9   "
                              "-4.544199297801071e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 8      2 10    "
                              "6.985326764982149e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8       2 2    "
                              "2.742839545884835e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8       2 3   "
                              "-1.007254829333754e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8       2 4   "
                              "-4.238545057657446e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8       2 5   "
                              "-1.959727484419270e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8       2 6   "
                              "-4.359343171640759e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8       2 7    "
                              "2.870595717144795e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8       2 8    "
                              "2.858295831875154e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8       2 9   "
                              "-2.245129471809810e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 8      2 10    "
                              "3.366992515594730e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 1   "
                              "-1.988457274011615e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 2    "
                              "5.774732962362260e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 3   "
                              "-2.118969104465797e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 4   "
                              "-8.922444358257370e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 5   "
                              "-4.107677408558212e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 6   "
                              "-9.167119239052863e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 7    "
                              "6.038952957782194e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 8    "
                              "6.012286414752545e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8       2 9   "
                              "-4.723547742121070e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 8      2 10    "
                              "7.077437967354231e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 1    "
                              "1.746843528431625e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 2   "
                              "-5.061222867872834e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 3    "
                              "1.862452156449607e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 4    "
                              "7.824021992367237e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 5    "
                              "3.658743194107578e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 6    "
                              "8.069562405116350e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 7   "
                              "-5.308229091527655e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 8   "
                              "-5.287106143354928e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8       2 9    "
                              "4.150843319352120e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 8      2 10   "
                              "-6.239469776425440e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 8       2 8    "
                              "3.330593692195403e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 8       2 3   "
                              "-2.459244785532941e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 8       2 4   "
                              "-1.012354241756131e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 8       2 5   "
                              "-5.178464362429600e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 8       2 6   "
                              "-1.065161737061373e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 8       2 7    "
                              "6.929917642364796e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 8       2 8    "
                              "6.936027991558963e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 8       2 9   "
                              "-5.401560295211441e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 8      2 10    "
                              "8.287079967044912e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 9       2 5   "
                              "-2.641232153392844e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 9       2 6    "
                              "2.941297609416028e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 9       2 7   "
                              "-1.191198100873518e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 9       2 8    "
                              "1.261371567272619e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 9       2 9    "
                              "1.284858463010085e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 0       2 9      2 10    "
                              "5.517456553623749e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 9       2 3    "
                              "1.040082899463006e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 9       2 4   "
                              "-3.180643948960411e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 9       2 5    "
                              "9.607676452642700e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 9       2 6   "
                              "-1.068556733649408e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 9       2 7    "
                              "4.300596499617235e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 9       2 8   "
                              "-4.544199297801071e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 9       2 9   "
                              "-4.650551535095815e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2       2 9      2 10   "
                              "-1.979305449693643e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9       2 2    "
                              "1.789233425378183e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9       2 3    "
                              "5.219917234123776e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9       2 4   "
                              "-1.593432574526110e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9       2 5    "
                              "4.806215534947030e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9       2 6   "
                              "-5.364012338293058e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9       2 7    "
                              "2.122274753797954e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9       2 8   "
                              "-2.245129471809810e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9       2 9   "
                              "-2.312903290591722e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 4       2 9      2 10   "
                              "-9.689036069479793e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9       2 2    "
                              "3.766287808970467e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9       2 3    "
                              "1.098827415807882e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9       2 4   "
                              "-3.354025486260517e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9       2 5    "
                              "1.011625124175042e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9       2 6   "
                              "-1.129153857177379e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9       2 7    "
                              "4.465101923264921e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9       2 8   "
                              "-4.723547742121070e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9       2 9   "
                              "-4.867332297090299e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6       2 9      2 10   "
                              "-2.037840888536811e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9       2 2   "
                              "-3.303255145773652e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9       2 3   "
                              "-9.637227239216609e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9       2 4    "
                              "2.942090964056192e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9       2 5   "
                              "-8.876029555580901e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9       2 6    "
                              "9.902405437257963e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9       2 7   "
                              "-3.924516205341246e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9       2 8    "
                              "4.150843319352119e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9       2 9    "
                              "4.273668554975398e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 8       2 9      2 10    "
                              "1.792916316377662e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 9       2 9       2 9   "
                              "-2.651444217584578e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 9       2 3    "
                              "1.228450458146215e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 9       2 4   "
                              "-3.787129376070253e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 9       2 5    "
                              "1.135705037198195e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 9       2 6   "
                              "-1.269762453613440e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 9       2 7    "
                              "5.083215609747401e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 9       2 8   "
                              "-5.401560295211441e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 9       2 9   "
                              "-5.513634585383528e+02 \n";
    integralFileTwoBodyFAD << "       1 6      1 10       2 9      2 10   "
                              "-2.350614110853461e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0      2 10       2 4   "
                              "-1.086434005727397e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0      2 10       2 6   "
                              "-4.487282447436252e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0      2 10       2 7    "
                              "2.147959907634374e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0      2 10       2 8   "
                              "-1.959826459762755e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 0      2 10       2 9    "
                              "5.517456553623750e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 0      2 10      2 10   "
                              "-2.160981958127234e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 2      2 10       2 4    "
                              "3.905523455266483e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2      2 10       2 5   "
                              "-1.480250734659060e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2      2 10       2 6    "
                              "1.607878057938481e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2      2 10       2 7   "
                              "-7.750490431838681e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 2      2 10       2 8    "
                              "6.985326764982152e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2      2 10       2 9   "
                              "-1.979305449693643e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 2      2 10      2 10    "
                              "7.737309656852620e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10       2 2   "
                              "-1.825769883946855e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10       2 3    "
                              "2.886987091531489e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10       2 4    "
                              "1.929757887422544e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10       2 5   "
                              "-7.551003996121390e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10       2 6    "
                              "7.891806030829301e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10       2 7   "
                              "-3.850153106984002e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10       2 8    "
                              "3.366992515594731e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10       2 9   "
                              "-9.689036069479793e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 4      2 10      2 10    "
                              "3.780849035335164e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10       2 2   "
                              "-3.842278465413999e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10       2 3    "
                              "6.069100412154365e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10       2 4    "
                              "4.059993790292226e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10       2 5   "
                              "-1.590427668219145e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10       2 6    "
                              "1.659954508437191e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10       2 7   "
                              "-8.101896372023505e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10       2 8    "
                              "7.077437967354226e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10       2 9   "
                              "-2.037840888536811e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 6      2 10      2 10    "
                              "7.951316818881516e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 7      2 10      2 10   "
                              "-1.280994572392286e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10       2 2    "
                              "3.372852130492381e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10       2 3   "
                              "-5.348072594774358e+00 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10       2 4   "
                              "-3.567575078941885e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10       2 5    "
                              "1.391957884620223e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10       2 6   "
                              "-1.459885354939897e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10       2 7    "
                              "7.114493677598676e+01 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10       2 8   "
                              "-6.239469776425442e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10       2 9    "
                              "1.792916316377662e+02 \n";
    integralFileTwoBodyFAD << "       1 6       1 8      2 10      2 10   "
                              "-6.997127383543939e+03 \n";
    integralFileTwoBodyFAD << "       1 6       1 9      2 10      2 10    "
                              "4.422719749063041e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10      2 10       2 4    "
                              "4.660740009072177e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10      2 10       2 5   "
                              "-1.757020148799409e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10      2 10       2 6    "
                              "1.918822020574652e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10      2 10       2 7   "
                              "-9.231586122277202e+00 \n";
    integralFileTwoBodyFAD << "       1 6      1 10      2 10       2 8    "
                              "8.287079967044912e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10      2 10       2 9   "
                              "-2.350614110853462e+01 \n";
    integralFileTwoBodyFAD << "       1 6      1 10      2 10      2 10    "
                              "9.211453258183600e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 0       2 1   "
                              "-2.895320936801981e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 0       2 2    "
                              "2.274405728380949e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 0       2 1   "
                              "-1.548304605493568e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 0       2 2    "
                              "1.193959336901426e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 0       2 0    "
                              "3.586983501749820e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 0       2 1   "
                              "-5.596348316474127e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 0       2 2    "
                              "4.223994946478371e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 0       2 3    "
                              "3.089803528482303e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 0       2 0   "
                              "-1.291208867215283e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 0       2 1    "
                              "1.863872830025212e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 0       2 2   "
                              "-1.393588908891451e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 0       2 3   "
                              "-9.752394052974573e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 0       2 4   "
                              "-1.703314529560679e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 0       2 0    "
                              "1.799256995632456e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 0       2 1   "
                              "-3.357053977900662e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 0       2 2    "
                              "2.581993335835812e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 0       2 3    "
                              "1.961594360186808e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 0       2 1    "
                              "1.503299029713919e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 1       2 0   "
                              "-2.895320936801981e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 1       2 2    "
                              "4.099039183720793e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 1       2 3    "
                              "3.870940551515010e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 1       2 0   "
                              "-1.548304605493568e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 1       2 1   "
                              "-1.052824373302266e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 1       2 2    "
                              "2.191766868201736e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 1       2 3    "
                              "2.033046895459608e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 1       2 4    "
                              "3.209946277326264e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 1       2 0   "
                              "-5.596348316474127e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 1       2 1   "
                              "-1.836616689658290e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 1       2 2    "
                              "7.918326825421234e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 1       2 3    "
                              "7.193974856486133e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 1       2 4    "
                              "1.092206969244391e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 1       2 5   "
                              "-1.468160434437459e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 1       2 6    "
                              "2.389380046566018e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 1       2 8    "
                              "1.081104317671464e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 0    "
                              "1.863872830025212e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 1    "
                              "3.274745130764540e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 2   "
                              "-2.636534013667273e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 3   "
                              "-2.373520397088738e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 4   "
                              "-3.524131616410955e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 5    "
                              "4.836641027457873e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 6   "
                              "-8.015107433658073e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 7    "
                              "2.803263633480148e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 8   "
                              "-3.593056323767811e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1       2 9    "
                              "1.147867459059400e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 1      2 10    "
                              "1.580396119094091e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 1       2 0   "
                              "-3.357053977900662e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 1       2 1   "
                              "-2.137894636141066e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 1       2 2    "
                              "4.750875450135381e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 1       2 3    "
                              "4.395271864661642e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 1       2 4    "
                              "6.781248699377402e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 1       2 6    "
                              "1.435419026291358e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 1       2 0    "
                              "1.503299029713919e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 1       2 2   "
                              "-2.127438457360215e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 2       2 0    "
                              "2.274405728380951e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 2       2 1    "
                              "4.099039183720793e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 2       2 3   "
                              "-4.954502197377936e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 2       2 4   "
                              "-6.199420936168711e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 2       2 0    "
                              "1.193959336901426e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 2       2 1    "
                              "2.191766868201736e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 2       2 2   "
                              "-3.149375302594184e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 2       2 3   "
                              "-2.648950877321754e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 2       2 4   "
                              "-3.262191804448531e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 2       2 7   "
                              "-1.125542857857265e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 0    "
                              "4.223994946478371e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 1    "
                              "7.918326825421234e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 2   "
                              "-8.209974031941123e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 3   "
                              "-9.565319590112397e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 4   "
                              "-1.156118632219219e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 5   "
                              "-3.349919430076870e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 7   "
                              "-4.080648406830742e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 8   "
                              "-3.119557424450925e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2       2 9   "
                              "-2.038515382809466e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 2      2 10    "
                              "2.084417658589885e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 0   "
                              "-1.393588908891451e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 1   "
                              "-2.636534013667273e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 2    "
                              "2.275154907966963e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 3    "
                              "3.184080623768809e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 4    "
                              "3.816651787518476e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 5    "
                              "1.246667732262860e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 7    "
                              "1.361780590931601e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 8    "
                              "1.043346691893512e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2       2 9    "
                              "6.804967209970212e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 2      2 10   "
                              "-6.942717033544527e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 0    "
                              "2.581993335835809e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 1    "
                              "4.750875450135381e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 2   "
                              "-6.583433655063602e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 3   "
                              "-5.740122322743413e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 4   "
                              "-7.050524343971142e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 5   "
                              "-1.756238495991521e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 7   "
                              "-2.448791263423903e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 8   "
                              "-1.861023223998564e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2       2 9   "
                              "-1.220413371137132e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 2      2 10    "
                              "1.251596804748730e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 2       2 1   "
                              "-2.127438457360214e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 2       2 3    "
                              "2.570381918532748e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 3       2 1    "
                              "3.870940551515007e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 3       2 2   "
                              "-4.954502197377937e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 3       2 4   "
                              "-6.004325429104954e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 3       2 5    "
                              "2.708858206606645e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 3       2 6    "
                              "4.534297069014698e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 1    "
                              "2.033046895459609e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 2   "
                              "-2.648950877321754e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 3    "
                              "4.752201350175850e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 4   "
                              "-3.209341015703792e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 5    "
                              "1.395360089645935e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 6    "
                              "2.396859159788276e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 7   "
                              "-4.816566679604800e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 8    "
                              "3.221609331654035e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 3       2 9   "
                              "-1.631812734709109e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 0    "
                              "3.089803528482288e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 1    "
                              "7.193974856486135e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 2   "
                              "-9.565319590112397e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 3    "
                              "2.171104921565691e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 4   "
                              "-1.158248083248140e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 5    "
                              "4.841163461759433e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 6    "
                              "8.459141036347975e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 7   "
                              "-1.757941117076523e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 8    "
                              "1.153836792718764e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3       2 9   "
                              "-5.943802001611488e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 3      2 10   "
                              "-3.325931811511966e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 0   "
                              "-9.752394052974402e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 1   "
                              "-2.373520397088737e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 2    "
                              "3.184080623768809e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 3   "
                              "-7.888213081336143e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 4    "
                              "3.854504161906462e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 5   "
                              "-1.583723886747646e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 6   "
                              "-2.783772333046887e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 7    "
                              "5.877121273357159e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 8   "
                              "-3.828897218869155e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3       2 9    "
                              "1.985579965391311e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 3      2 10    "
                              "1.096931156114034e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 0    "
                              "1.961594360186808e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 1    "
                              "4.395271864661641e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 2   "
                              "-5.740122322743413e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 3    "
                              "1.066770086804327e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 4   "
                              "-6.952975158208280e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 5    "
                              "3.014409335018217e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 6    "
                              "5.147142479302507e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 7   "
                              "-1.043490069029160e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 8    "
                              "6.990737527840296e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3       2 9   "
                              "-3.536361341968045e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 3      2 10   "
                              "-2.024942526576513e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 3       2 2    "
                              "2.570381918532748e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 3       2 4    "
                              "3.113584990697548e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 4       2 2   "
                              "-6.199420936168713e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 4       2 3   "
                              "-6.004325429104954e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 4       2 4   "
                              "-1.091113269607886e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 4       2 5    "
                              "5.195899699454706e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 4       2 6    "
                              "2.085247865990193e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 4       2 7    "
                              "2.426654210686757e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 4       2 8    "
                              "2.492485166653134e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 4      2 10   "
                              "-1.152958553230018e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 1    "
                              "3.209946277326271e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 2   "
                              "-3.262191804448531e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 3   "
                              "-3.209341015703792e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 4   "
                              "-5.715926739403925e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 5    "
                              "2.779598137705229e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 6    "
                              "1.106400789464097e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 7    "
                              "1.304669400881715e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 8    "
                              "1.330679042005104e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4       2 9    "
                              "5.002264664759648e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 4      2 10   "
                              "-6.142198838335279e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 1    "
                              "1.092206969244393e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 2   "
                              "-1.156118632219220e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 3   "
                              "-1.158248083248140e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 4   "
                              "-2.014944216893012e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 5    "
                              "1.003491365462847e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 6    "
                              "3.952745899522264e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 7    "
                              "4.781871010178237e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 8    "
                              "4.828240276279019e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4       2 9    "
                              "1.815838556879361e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 4      2 10   "
                              "-2.208292144964530e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 0   "
                              "-1.703314529560679e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 1   "
                              "-3.524131616410961e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 2    "
                              "3.816651787518476e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 3    "
                              "3.854504161906462e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 4    "
                              "6.635900382588818e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 5   "
                              "-3.339703421393776e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 6   "
                              "-1.309286385919521e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 7   "
                              "-1.603744021364651e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 8   "
                              "-1.612013728478309e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4       2 9   "
                              "-6.060278108033321e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 4      2 10    "
                              "7.336705740688581e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 1    "
                              "6.781248699377402e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 2   "
                              "-7.050524343971142e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 3   "
                              "-6.952975158208280e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 4   "
                              "-1.234468756169917e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 5    "
                              "6.019892800671735e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 6    "
                              "2.391210127227888e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 7    "
                              "2.846648138015005e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 8    "
                              "2.898027478814581e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4       2 9    "
                              "1.086284764807786e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 4      2 10   "
                              "-1.330019087231059e+01 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 4       2 3    "
                              "3.113584990697549e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 4       2 5   "
                              "-2.695151589885668e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 4       2 6   "
                              "-1.072252078131341e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 5       2 3    "
                              "2.708858206606645e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 5       2 4    "
                              "5.195899699454706e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 5       2 5    "
                              "4.662981931487190e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 5       2 6   "
                              "-6.377433945415034e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 5       2 7   "
                              "-1.194218870577145e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 5       2 8    "
                              "1.316093853223485e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 5       2 9   "
                              "-2.799012360995763e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 5       2 3    "
                              "1.395360089645935e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 5       2 4    "
                              "2.779598137705229e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 5       2 5    "
                              "2.509255609007438e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 5       2 6   "
                              "-3.405310144163745e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 5       2 7   "
                              "-6.434589495140034e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 5       2 8    "
                              "6.840243717938973e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 5       2 9   "
                              "-1.507482553964324e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 5      2 10    "
                              "2.327909692629904e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 1   "
                              "-1.468160434437445e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 2   "
                              "-3.349919430076785e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 3    "
                              "4.841163461759429e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 4    "
                              "1.003491365462848e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 5    "
                              "9.124093220817822e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 6   "
                              "-1.226713304451232e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 7   "
                              "-2.333108475191311e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 8    "
                              "2.311902370230600e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5       2 9   "
                              "-5.478109977465257e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 5      2 10    "
                              "8.558681798928369e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 1    "
                              "4.836641027457930e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 2    "
                              "1.246667732262852e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 3   "
                              "-1.583723886747648e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 4   "
                              "-3.339703421393776e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 5   "
                              "-3.046234584089571e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 6    "
                              "4.078700818095911e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 7    "
                              "7.775807930212953e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 8   "
                              "-7.431598220832090e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5       2 9    "
                              "1.828050525580669e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 5      2 10   "
                              "-2.873274308630607e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5       2 2   "
                              "-1.756238495991539e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5       2 3    "
                              "3.014409335018217e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5       2 4    "
                              "6.019892800671736e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5       2 5    "
                              "5.441386129562741e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5       2 6   "
                              "-7.372631609575187e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5       2 7   "
                              "-1.391514447957065e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5       2 8    "
                              "1.439431409543724e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5       2 9   "
                              "-3.265023418343500e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 5      2 10    "
                              "5.074491052839496e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 5       2 4   "
                              "-2.695151589885669e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 5       2 5   "
                              "-2.433439865912350e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 5       2 6    "
                              "3.301975296170514e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 6       2 3    "
                              "4.534297069014702e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 6       2 4    "
                              "2.085247865990194e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 6       2 5   "
                              "-6.377433945415034e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 6       2 6   "
                              "-8.689957917015737e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 6       2 7    "
                              "2.744790145473134e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 6       2 8    "
                              "2.645144412229411e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 6       2 9    "
                              "3.127256675103661e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 6      2 10   "
                              "-4.755186837855244e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 6       2 3    "
                              "2.396859159788276e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 6       2 4    "
                              "1.106400789464097e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 6       2 5   "
                              "-3.405310144163745e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 6       2 6   "
                              "-4.625384216280532e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 6       2 7    "
                              "1.469152080985388e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 6       2 8    "
                              "1.404713045722655e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 6       2 9    "
                              "1.679625509592282e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 6      2 10   "
                              "-2.527832710654950e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6       2 1    "
                              "2.389380046566017e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6       2 3    "
                              "8.459141036347968e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6       2 4    "
                              "3.952745899522264e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6       2 5   "
                              "-1.226713304451232e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6       2 6   "
                              "-1.660072126397794e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6       2 7    "
                              "5.302471945312025e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6       2 8    "
                              "5.008491166851751e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6       2 9    "
                              "6.108645948415147e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 6      2 10   "
                              "-9.048897347514762e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6       2 1   "
                              "-8.015107433657960e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6       2 3   "
                              "-2.783772333046887e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6       2 4   "
                              "-1.309286385919521e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6       2 5    "
                              "4.078700818095911e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6       2 6    "
                              "5.510575254689214e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6       2 7   "
                              "-1.764295338504232e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6       2 8   "
                              "-1.656828356052032e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6       2 9   "
                              "-2.040240376923094e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 6      2 10    "
                              "2.999865254453595e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6       2 1    "
                              "1.435419026291357e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6       2 3    "
                              "5.147142479302509e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6       2 4    "
                              "2.391210127227888e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6       2 5   "
                              "-7.372631609575187e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6       2 6   "
                              "-1.000866563941850e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6       2 7    "
                              "3.180067623445261e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6       2 8    "
                              "3.029150003666390e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6       2 9    "
                              "3.647176483226019e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 6      2 10   "
                              "-5.464381642888385e+01 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 6       2 4   "
                              "-1.072252078131340e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 6       2 5    "
                              "3.301975296170513e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 6       2 6    "
                              "4.485429484788905e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 6       2 7   "
                              "-1.423405578005329e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 6       2 8   "
                              "-1.357942837917770e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 7       2 4    "
                              "2.426654210686757e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 7       2 5   "
                              "-1.194218870577144e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 7       2 6    "
                              "2.744790145473134e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 7       2 7    "
                              "1.075387200487358e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 7       2 8   "
                              "-1.715197084999769e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 7       2 9   "
                              "-1.256027146080468e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 7      2 10    "
                              "2.279602297199011e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7       2 2   "
                              "-1.125542857857275e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7       2 3   "
                              "-4.816566679604801e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7       2 4    "
                              "1.304669400881717e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7       2 5   "
                              "-6.434589495140034e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7       2 6    "
                              "1.469152080985388e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7       2 7    "
                              "5.771349794646540e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7       2 8   "
                              "-9.149074069379840e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7       2 9   "
                              "-6.742517741361078e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 7      2 10    "
                              "1.218472311167989e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7       2 2   "
                              "-4.080648406830771e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7       2 3   "
                              "-1.757941117076527e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7       2 4    "
                              "4.781871010178237e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7       2 5   "
                              "-2.333108475191311e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7       2 6    "
                              "5.302471945312025e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7       2 7    "
                              "2.090162750628616e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7       2 8   "
                              "-3.286902251655388e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7       2 9   "
                              "-2.428842699671646e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 7      2 10    "
                              "4.398595092979115e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 1    "
                              "2.803263633480113e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 2    "
                              "1.361780590931609e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 3    "
                              "5.877121273357171e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 4   "
                              "-1.603744021364654e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 5    "
                              "7.775807930212953e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 6   "
                              "-1.764295338504231e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 7   "
                              "-6.965410046168022e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 8    "
                              "1.091342497669830e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7       2 9    "
                              "8.069740539693474e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 7      2 10   "
                              "-1.464027450350100e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7       2 2   "
                              "-2.448791263423904e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7       2 3   "
                              "-1.043490069029159e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7       2 4    "
                              "2.846648138015014e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7       2 5   "
                              "-1.391514447957065e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7       2 6    "
                              "3.180067623445261e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7       2 7    "
                              "1.250162959060689e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7       2 8   "
                              "-1.978191421576763e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7       2 9   "
                              "-1.455209600846862e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 7      2 10    "
                              "2.640636332093960e+01 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 7       2 6   "
                              "-1.423405578005329e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 7       2 7   "
                              "-5.592939369228017e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 8       2 4    "
                              "2.492485166653136e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 8       2 5    "
                              "1.316093853223485e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 8       2 6    "
                              "2.645144412229411e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 8       2 7   "
                              "-1.715197084999769e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 8       2 8   "
                              "-1.718651925716282e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 8       2 9    "
                              "1.335862154989680e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 8      2 10   "
                              "-2.064903844978515e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 8       2 3    "
                              "3.221609331654021e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 8       2 4    "
                              "1.330679042005102e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 8       2 5    "
                              "6.840243717938989e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 8       2 6    "
                              "1.404713045722655e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 8       2 7   "
                              "-9.149074069379840e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 8       2 8   "
                              "-9.145112929589289e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 8       2 9    "
                              "7.137817176133413e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 8      2 10   "
                              "-1.095365014311530e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 1    "
                              "1.081104317671464e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 2   "
                              "-3.119557424450922e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 3    "
                              "1.153836792718764e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 4    "
                              "4.828240276279025e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 5    "
                              "2.311902370230600e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 6    "
                              "5.008491166851751e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 7   "
                              "-3.286902251655388e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 8   "
                              "-3.276419930057813e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8       2 9    "
                              "2.568935552594200e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 8      2 10   "
                              "-3.881011585810597e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 1   "
                              "-3.593056323767839e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 2    "
                              "1.043346691893501e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 3   "
                              "-3.828897218869161e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 4   "
                              "-1.612013728478307e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 5   "
                              "-7.431598220832090e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 6   "
                              "-1.656828356052032e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 7    "
                              "1.091342497669830e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 8    "
                              "1.086541346241882e+04 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8       2 9   "
                              "-8.536183623767128e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 8      2 10    "
                              "1.279349171509797e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 8       2 8       2 8   "
                              "-1.235972933403576e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8       2 2   "
                              "-1.861023223998550e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8       2 3    "
                              "6.990737527840335e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8       2 4    "
                              "2.898027478814583e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8       2 5    "
                              "1.439431409543724e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8       2 6    "
                              "3.029150003666390e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8       2 7   "
                              "-1.978191421576763e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8       2 8   "
                              "-1.976566944719781e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8       2 9    "
                              "1.543651779489834e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 8      2 10   "
                              "-2.352073556074062e+02 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 8       2 6   "
                              "-1.357942837917770e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 8       2 8    "
                              "8.858111689187398e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 8      2 10    "
                              "1.054592973375672e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 9       2 5   "
                              "-2.799012360995763e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 9       2 6    "
                              "3.127256675103663e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 9       2 7   "
                              "-1.256027146080468e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 9       2 8    "
                              "1.335862154989680e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 9       2 9    "
                              "1.360631752145197e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 1       2 9      2 10    "
                              "5.826741370141415e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 9       2 3   "
                              "-1.631812734709112e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 9       2 4    "
                              "5.002264664759646e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 9       2 5   "
                              "-1.507482553964327e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 9       2 6    "
                              "1.679625509592283e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 9       2 7   "
                              "-6.742517741361075e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 9       2 8    "
                              "7.137817176133413e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 9       2 9    "
                              "7.301264573074578e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3       2 9      2 10    "
                              "3.106728481330243e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9       2 2   "
                              "-2.038515382809466e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9       2 3   "
                              "-5.943802001611488e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9       2 4    "
                              "1.815838556879361e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9       2 5   "
                              "-5.478109977465257e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9       2 6    "
                              "6.108645948415139e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9       2 7   "
                              "-2.428842699671647e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9       2 8    "
                              "2.568935552594200e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9       2 9    "
                              "2.641087236144967e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5       2 9      2 10    "
                              "1.111693254863085e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 1    "
                              "1.147867459059422e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 2    "
                              "6.804967209970195e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 3    "
                              "1.985579965391311e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 4   "
                              "-6.060278108033338e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 5    "
                              "1.828050525580664e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 6   "
                              "-2.040240376923094e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 7    "
                              "8.069740539693488e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 8   "
                              "-8.536183623767128e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9       2 9   "
                              "-8.795696732602295e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7       2 9      2 10   "
                              "-3.683071836158664e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9       2 2   "
                              "-1.220413371137132e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9       2 3   "
                              "-3.536361341968047e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9       2 4    "
                              "1.086284764807787e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9       2 5   "
                              "-3.265023418343506e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9       2 6    "
                              "3.647176483226018e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9       2 7   "
                              "-1.455209600846862e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9       2 8    "
                              "1.543651779489834e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9       2 9    "
                              "1.580454111766437e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 9       2 9      2 10    "
                              "6.700652879212963e+01 \n";
    integralFileTwoBodyFAD << "       1 7      1 10       2 9       2 9   "
                              "-7.069141743283788e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1      2 10       2 4   "
                              "-1.152958553230019e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1      2 10       2 6   "
                              "-4.755186837855244e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1      2 10       2 7    "
                              "2.279602297199010e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1      2 10       2 8   "
                              "-2.064903844978516e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 1      2 10       2 9    "
                              "5.826741370141416e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 1      2 10      2 10   "
                              "-2.286064961767198e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3      2 10       2 4   "
                              "-6.142198838335279e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3      2 10       2 5    "
                              "2.327909692629904e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 3      2 10       2 6   "
                              "-2.527832710654951e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3      2 10       2 7    "
                              "1.218472311167989e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3      2 10       2 8   "
                              "-1.095365014311529e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 3      2 10       2 9    "
                              "3.106728481330242e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 3      2 10      2 10   "
                              "-1.215325755672726e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10       2 2    "
                              "2.084417658589888e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10       2 3   "
                              "-3.325931811511966e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10       2 4   "
                              "-2.208292144964530e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10       2 5    "
                              "8.558681798928353e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10       2 6   "
                              "-9.048897347514762e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10       2 7    "
                              "4.398595092979114e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10       2 8   "
                              "-3.881011585810597e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10       2 9    "
                              "1.111693254863085e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 5      2 10      2 10   "
                              "-4.340608353924198e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 6      2 10      2 10   "
                              "-1.280994572393041e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 1    "
                              "1.580396119094087e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 2   "
                              "-6.942717033544541e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 3    "
                              "1.096931156114040e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 4    "
                              "7.336705740688573e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 5   "
                              "-2.873274308630601e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 6    "
                              "2.999865254453595e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 7   "
                              "-1.464027450350100e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 8    "
                              "1.279349171509797e+03 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10       2 9   "
                              "-3.683071836158663e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 7      2 10      2 10    "
                              "1.437059650309869e+04 \n";
    integralFileTwoBodyFAD << "       1 7       1 8      2 10      2 10   "
                              "-1.653820579459914e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10       2 2    "
                              "1.251596804748731e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10       2 3   "
                              "-2.024942526576514e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10       2 4   "
                              "-1.330019087231059e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10       2 5    "
                              "5.074491052839496e+00 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10       2 6   "
                              "-5.464381642888385e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10       2 7    "
                              "2.640636332093960e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10       2 8   "
                              "-2.352073556074063e+02 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10       2 9    "
                              "6.700652879212963e+01 \n";
    integralFileTwoBodyFAD << "       1 7       1 9      2 10      2 10   "
                              "-2.621992222462075e+03 \n";
    integralFileTwoBodyFAD << "       1 7      1 10      2 10       2 8    "
                              "1.054592973375673e+00 \n";
    integralFileTwoBodyFAD << "       1 7      1 10      2 10      2 10    "
                              "1.175238565629099e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 0       2 1    "
                              "9.898910128722843e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 0       2 1   "
                              "-3.151543766714967e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 0       2 2    "
                              "2.481706819616808e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 0       2 1   "
                              "-1.330708936734694e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 0       2 2    "
                              "1.042637926696740e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 0       2 0    "
                              "6.111356830875017e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 0       2 1   "
                              "-9.055431617871011e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 0       2 2    "
                              "6.792768343223963e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 0       2 3    "
                              "4.834720798096543e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 0       2 0   "
                              "-8.208710741855775e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 0       2 1    "
                              "1.166483838714265e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 0       2 2   "
                              "-8.704110725180921e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 0       2 3   "
                              "-6.023150268082939e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 0       2 4   "
                              "-1.052987922447612e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 0       2 1   "
                              "-1.609511939947187e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 0       2 0    "
                              "1.992058670279731e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 0       2 1   "
                              "-3.600904551175626e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 0       2 2    "
                              "2.761058304293002e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 0       2 3    "
                              "2.082722879847968e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 1       2 0    "
                              "9.898910128722845e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 1       2 2   "
                              "-1.400968750064177e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 1       2 3   "
                              "-1.326051652580064e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 1       2 0   "
                              "-3.151543766714967e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 1       2 2    "
                              "4.461175351248307e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 1       2 3    "
                              "4.222429427656756e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 1       2 0   "
                              "-1.330708936734694e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 1       2 1   "
                              "-1.259122793097391e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 1       2 2    "
                              "1.884671390058141e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 1       2 3    "
                              "1.775379814730944e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 1       2 4    "
                              "2.910704981605829e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 0   "
                              "-9.055431617871011e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 1   "
                              "-2.067448858236699e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 2    "
                              "1.281057154462392e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 3    "
                              "1.156926370480235e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 4    "
                              "1.732519626425011e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 5   "
                              "-2.361289543442362e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 6    "
                              "3.882856635079146e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 7   "
                              "-1.359032754209603e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 1       2 8    "
                              "1.746843528431626e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 0    "
                              "1.166483838714265e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 1    "
                              "1.673092278198236e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 2   "
                              "-1.649939878900873e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 3   "
                              "-1.482448735576118e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 4   "
                              "-2.188630503399878e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 5    "
                              "3.016277859113235e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 6   "
                              "-5.026060072298471e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 7    "
                              "1.756686374900237e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 1       2 8   "
                              "-2.247776729083665e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 1       2 0   "
                              "-1.609511939947187e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 1       2 2    "
                              "2.277692492269826e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 1       2 0   "
                              "-3.600904551175626e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 1       2 1   "
                              "-2.110613471131270e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 1       2 2    "
                              "5.095763963044434e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 1       2 3    "
                              "4.700424631374047e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 1       2 4    "
                              "7.228580193993118e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 1       2 6    "
                              "1.539879421310958e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 2       2 1   "
                              "-1.400968750064177e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 2       2 3    "
                              "1.692724307676117e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 2       2 4    "
                              "2.122050539075623e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 2       2 0    "
                              "2.481706819616808e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 2       2 1    "
                              "4.461175351248307e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 2       2 3   "
                              "-5.391394663410122e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 2       2 4   "
                              "-6.758362232961090e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 2       2 0    "
                              "1.042637926696741e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 2       2 1    "
                              "1.884671390058140e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 2       2 2   "
                              "-3.279620017120516e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 2       2 3   "
                              "-2.278966218560902e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 2       2 4   "
                              "-2.846256946553360e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 0    "
                              "6.792768343223963e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 1    "
                              "1.281057154462392e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 2   "
                              "-1.182443949774317e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 3   "
                              "-1.547259245000407e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 4   "
                              "-1.860006153611544e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 5   "
                              "-5.821056516133534e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 7   "
                              "-6.610698557937487e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 8   "
                              "-5.061222867872859e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2       2 9   "
                              "-3.303255145773659e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 2      2 10    "
                              "3.372852130492381e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 0   "
                              "-8.704110725180924e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 1   "
                              "-1.649939878900873e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 2    "
                              "1.362913529394083e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 3    "
                              "1.992461346737492e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 4    "
                              "2.384049076928480e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 5    "
                              "7.998397159617324e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 7    "
                              "8.527275601818486e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 8    "
                              "6.536053782709580e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2       2 9    "
                              "4.261223997535487e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 2      2 10   "
                              "-4.345238508775916e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 2       2 1    "
                              "2.277692492269826e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 2       2 3   "
                              "-2.751898147320738e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 0    "
                              "2.761058304293000e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 1    "
                              "5.095763963044434e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 2   "
                              "-6.769114163844153e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 3   "
                              "-6.156582684296372e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 4   "
                              "-7.542149448359862e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 5   "
                              "-1.934359721557124e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 7   "
                              "-2.626784858720455e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 8   "
                              "-1.998170202978411e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2       2 9   "
                              "-1.309575084929296e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 2      2 10    "
                              "1.342331337795866e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 3       2 1   "
                              "-1.326051652580065e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 3       2 2    "
                              "1.692724307676117e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 3       2 4    "
                              "2.051003749055948e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 3       2 6   "
                              "-1.537535523632358e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 3       2 1    "
                              "4.222429427656756e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 3       2 2   "
                              "-5.391394663410124e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 3       2 4   "
                              "-6.533408977094894e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 3       2 5    "
                              "2.964459475455942e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 3       2 6    "
                              "4.920829119399282e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 1    "
                              "1.775379814730944e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 2   "
                              "-2.278966218560902e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 3    "
                              "3.263834189772104e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 4   "
                              "-2.762467533221910e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 5    "
                              "1.234450371072647e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 6    "
                              "2.104573333345820e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 7   "
                              "-4.113180474383497e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 8    "
                              "2.785479656651536e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 3       2 9   "
                              "-1.395412567391564e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 0    "
                              "4.834720798096551e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 1    "
                              "1.156926370480236e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 2   "
                              "-1.547259245000407e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 3    "
                              "3.721939784918222e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 4   "
                              "-1.873227694276366e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 5    "
                              "7.741903802694505e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 6    "
                              "1.358448459682019e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 7   "
                              "-2.851818295367525e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 8    "
                              "1.862452156449607e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3       2 9   "
                              "-9.637227239216637e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 3      2 10   "
                              "-5.348072594774349e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 0   "
                              "-6.023150268082939e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 1   "
                              "-1.482448735576117e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 2    "
                              "1.992461346737492e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 3   "
                              "-5.024206442267462e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 4    "
                              "2.411830621959173e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 5   "
                              "-9.874247336925644e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 6   "
                              "-1.737244197191931e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 7    "
                              "3.680775561044462e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 8   "
                              "-2.394602011086547e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3       2 9    "
                              "1.243370880031011e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 3      2 10    "
                              "6.849785581267615e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 3       2 2   "
                              "-2.751898147320738e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 3       2 4   "
                              "-3.333024663302277e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 0    "
                              "2.082722879847989e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 1    "
                              "4.700424631374044e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 2   "
                              "-6.156582684296372e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 3    "
                              "1.185843945031393e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 4   "
                              "-7.456977486817235e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 5    "
                              "3.214129084211645e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 6    "
                              "5.506918924840453e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 7   "
                              "-1.121168330839419e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 8    "
                              "7.486746703671655e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3       2 9   "
                              "-3.798183267541381e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 3      2 10   "
                              "-2.166530796206881e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 4       2 2    "
                              "2.122050539075623e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 4       2 3    "
                              "2.051003749055948e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 4       2 4    "
                              "3.737094920055658e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 4       2 5   "
                              "-1.773536430507334e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 4       2 6   "
                              "-7.119597352753558e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4       2 2   "
                              "-6.758362232961091e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4       2 3   "
                              "-6.533408977094894e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4       2 4   "
                              "-1.190083045446680e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4       2 5    "
                              "5.651696464201537e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4       2 6    "
                              "2.269296754460614e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4       2 7    "
                              "2.646964005816338e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4       2 8    "
                              "2.721165044127955e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4       2 9    "
                              "1.016320947901380e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 4      2 10   "
                              "-1.255070378117837e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 1    "
                              "2.910704981605822e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 2   "
                              "-2.846256946553359e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 3   "
                              "-2.762467533221910e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 4   "
                              "-5.006662557819731e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 5    "
                              "2.392495473378872e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 6    "
                              "9.601837804309102e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 7    "
                              "1.106297935906884e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 8    "
                              "1.137894864898769e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4       2 9    "
                              "4.282577880577525e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 4      2 10   "
                              "-5.302150443587363e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 1    "
                              "1.732519626425011e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 2   "
                              "-1.860006153611544e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 3   "
                              "-1.873227694276366e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 4   "
                              "-3.236647582216826e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 5    "
                              "1.623031378217881e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 6    "
                              "6.373507062916849e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 7    "
                              "7.771610193093110e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 8    "
                              "7.824021992367236e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4       2 9    "
                              "2.942090964056198e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 4      2 10   "
                              "-3.567575078941882e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 0   "
                              "-1.052987922447670e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 1   "
                              "-2.188630503399878e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 2    "
                              "2.384049076928480e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 3    "
                              "2.411830621959174e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 4    "
                              "4.142928059430279e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 5   "
                              "-2.089704987287585e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 6   "
                              "-8.183890475657564e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 7   "
                              "-1.005349395803509e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 8   "
                              "-1.009548657688206e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4       2 9   "
                              "-3.794591879474167e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 4      2 10    "
                              "4.589079566073323e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 4       2 3   "
                              "-3.333024663302277e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 4       2 5    "
                              "2.886762933707565e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 4       2 6    "
                              "1.142896303026901e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 1    "
                              "7.228580193993140e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 2   "
                              "-7.542149448359862e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 3   "
                              "-7.456977486817235e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 4   "
                              "-1.319562657339695e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 5    "
                              "6.456887602996154e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 6    "
                              "2.561210028161616e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 7    "
                              "3.057786407931123e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 8    "
                              "3.108631840358291e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4       2 9    "
                              "1.165758245058048e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 4      2 10   "
                              "-1.425656523948272e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 5       2 4   "
                              "-1.773536430507334e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 5       2 5   "
                              "-1.591035244375781e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 5       2 6    "
                              "2.177318971603239e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 5       2 7    "
                              "4.054464594139068e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 5       2 3    "
                              "2.964459475455946e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 5       2 4    "
                              "5.651696464201537e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 5       2 5    "
                              "5.069600097409608e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 5       2 6   "
                              "-6.938272899615849e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 5       2 7   "
                              "-1.296239042406352e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 5       2 8    "
                              "1.417263682290600e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 5       2 9   "
                              "-3.039926879078881e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 5       2 3    "
                              "1.234450371072647e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 5       2 4    "
                              "2.392495473378872e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 5       2 5    "
                              "2.147640920859369e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 5       2 6   "
                              "-2.935943264637376e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 5       2 7   "
                              "-5.527434010713307e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 5       2 8    "
                              "6.237387858775524e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 5       2 9   "
                              "-1.291638838655700e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 5      2 10    "
                              "1.971049907245196e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 1   "
                              "-2.361289543442359e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 2   "
                              "-5.821056516133477e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 3    "
                              "7.741903802694516e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 4    "
                              "1.623031378217881e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 5    "
                              "1.478768439004187e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 6   "
                              "-1.982828794678772e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 7   "
                              "-3.777338787717429e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 8    "
                              "3.658743194107578e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5       2 9   "
                              "-8.876029555580901e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 5      2 10    "
                              "1.391957884620224e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 1    "
                              "3.016277859113150e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 2    "
                              "7.998397159617210e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 3   "
                              "-9.874247336925630e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 4   "
                              "-2.089704987287585e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 5   "
                              "-1.907389037980915e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 6    "
                              "2.551583884037253e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 7    "
                              "4.866493758373039e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 8   "
                              "-4.610877925549496e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5       2 9    "
                              "1.144451917307040e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 5      2 10   "
                              "-1.801466221487420e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 5       2 4    "
                              "2.886762933707565e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 5       2 5    "
                              "2.615391789411956e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 5       2 6   "
                              "-3.532842787495805e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5       2 2   "
                              "-1.934359721557124e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5       2 3    "
                              "3.214129084211645e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5       2 4    "
                              "6.456887602996154e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5       2 5    "
                              "5.842113970981338e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5       2 6   "
                              "-7.905409653530525e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5       2 7   "
                              "-1.493836520386780e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5       2 8    "
                              "1.533577798371569e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5       2 9   "
                              "-3.505700739697203e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 5      2 10    "
                              "5.454186991059930e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 6       2 3   "
                              "-1.537535523632358e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 6       2 4   "
                              "-7.119597352753558e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 6       2 5    "
                              "2.177318971603239e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 6       2 6    "
                              "2.968128860639843e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 6       2 7   "
                              "-9.359069647205237e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 6       2 8   "
                              "-9.005928584056248e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 6       2 9   "
                              "-1.070248767816951e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 6      2 10    "
                              "1.622634053486752e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 6       2 3    "
                              "4.920829119399281e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 6       2 4    "
                              "2.269296754460614e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 6       2 5   "
                              "-6.938272899615849e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 6       2 6   "
                              "-9.457325127999052e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 6       2 7    "
                              "2.984426762471246e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 6       2 8    "
                              "2.874383580633811e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 6       2 9    "
                              "3.403319184508235e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 6      2 10   "
                              "-5.173506711398748e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 6       2 3    "
                              "2.104573333345819e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 6       2 4    "
                              "9.601837804309102e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 6       2 5   "
                              "-2.935943264637376e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 6       2 6   "
                              "-3.999055142196679e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 6       2 7    "
                              "1.265206707577585e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 6       2 8    "
                              "1.222199276527281e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 6       2 9    "
                              "1.436182314881394e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 6      2 10   "
                              "-2.190682304157162e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6       2 1    "
                              "3.882856635079160e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6       2 3    "
                              "1.358448459682021e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6       2 4    "
                              "6.373507062916848e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6       2 5   "
                              "-1.982828794678772e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6       2 6   "
                              "-2.680443857574112e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6       2 7    "
                              "8.574994994802636e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6       2 8    "
                              "8.069562405116351e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6       2 9    "
                              "9.902405437257957e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 6      2 10   "
                              "-1.459885354939898e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6       2 1   "
                              "-5.026060072298440e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6       2 3   "
                              "-1.737244197191933e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6       2 4   "
                              "-8.183890475657564e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6       2 5    "
                              "2.551583884037253e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6       2 6    "
                              "3.446139044267455e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6       2 7   "
                              "-1.103868929624886e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6       2 8   "
                              "-1.035250899127595e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6       2 9   "
                              "-1.277659427864183e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 6      2 10    "
                              "1.875446118161310e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 6       2 4    "
                              "1.142896303026901e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 6       2 5   "
                              "-3.532842787495805e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 6       2 6   "
                              "-4.789834520304598e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 6       2 7    "
                              "1.525420447878652e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 6       2 8    "
                              "1.448191594496957e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6       2 1    "
                              "1.539879421310971e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6       2 3    "
                              "5.506918924840454e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6       2 4    "
                              "2.561210028161616e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6       2 5   "
                              "-7.905409653530525e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6       2 6   "
                              "-1.072637165170560e+03 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6       2 7    "
                              "3.411000185590887e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6       2 8    "
                              "3.244412860870672e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6       2 9    "
                              "3.915235751316374e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 6      2 10   "
                              "-5.854545589589042e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 7       2 5    "
                              "4.054464594139068e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 7       2 6   "
                              "-9.359069647205239e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 7       2 7   "
                              "-3.666759324278463e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 7       2 8    "
                              "5.848176844070158e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 7       2 9    "
                              "4.253366829644000e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 7       2 4    "
                              "2.646964005816338e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 7       2 5   "
                              "-1.296239042406352e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 7       2 6    "
                              "2.984426762471246e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 7       2 7    "
                              "1.169132270796727e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 7       2 8   "
                              "-1.865051437263648e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 7       2 9   "
                              "-1.363257048405140e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 7      2 10    "
                              "2.480996473328201e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 7       2 3   "
                              "-4.113180474383496e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 7       2 4    "
                              "1.106297935906884e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 7       2 5   "
                              "-5.527434010713307e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 7       2 6    "
                              "1.265206707577585e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 7       2 7    "
                              "4.956485384580502e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 7       2 8   "
                              "-7.908438434585076e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 7       2 9   "
                              "-5.825026979372014e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 7      2 10    "
                              "1.048492818872039e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 1   "
                              "-1.359032754209642e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 2   "
                              "-6.610698557937501e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 3   "
                              "-2.851818295367525e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 4    "
                              "7.771610193093110e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 5   "
                              "-3.777338787717429e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 6    "
                              "8.574994994802636e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 7    "
                              "3.383542788760048e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 8   "
                              "-5.308229091527654e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7       2 9   "
                              "-3.924516205341248e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 7      2 10    "
                              "7.114493677598679e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 1    "
                              "1.756686374900181e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 2    "
                              "8.527275601818427e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 3    "
                              "3.680775561044468e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 4   "
                              "-1.005349395803509e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 5    "
                              "4.866493758373039e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 6   "
                              "-1.103868929624886e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 7   "
                              "-4.359545575648399e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 8    "
                              "6.825000713784973e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7       2 9    "
                              "5.046870394011078e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 7      2 10   "
                              "-9.161076408066575e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 7       2 6    "
                              "1.525420447878651e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 7       2 7    "
                              "6.003160326643405e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7       2 2   "
                              "-2.626784858720455e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7       2 3   "
                              "-1.121168330839418e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7       2 4    "
                              "3.057786407931123e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7       2 5   "
                              "-1.493836520386780e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7       2 6    "
                              "3.411000185590887e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7       2 7    "
                              "1.341556434415410e+03 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7       2 8   "
                              "-2.120590029749823e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7       2 9   "
                              "-1.560969804035066e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 7      2 10    "
                              "2.832041378888437e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 8       2 6   "
                              "-9.005928584056251e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 8       2 7    "
                              "5.848176844070158e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 8       2 8    "
                              "5.862926376735704e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 8       2 9   "
                              "-4.552211551146172e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 8      2 10    "
                              "7.001489898205731e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 8       2 4    "
                              "2.721165044127956e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 8       2 5    "
                              "1.417263682290603e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 8       2 6    "
                              "2.874383580633811e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 8       2 7   "
                              "-1.865051437263648e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 8       2 8   "
                              "-1.869279983043623e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 8       2 9    "
                              "1.452317475658280e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 8      2 10   "
                              "-2.238895790085911e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 8       2 3    "
                              "2.785479656651535e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 8       2 4    "
                              "1.137894864898769e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 8       2 5    "
                              "6.237387858775537e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 8       2 6    "
                              "1.222199276527282e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 8       2 7   "
                              "-7.908438434585076e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 8       2 8   "
                              "-7.921312454946876e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 8       2 9    "
                              "6.161858483127648e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 8      2 10   "
                              "-9.591808561875580e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 1    "
                              "1.746843528431625e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 2   "
                              "-5.061222867872834e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 3    "
                              "1.862452156449607e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 4    "
                              "7.824021992367237e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 5    "
                              "3.658743194107578e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 6    "
                              "8.069562405116350e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 7   "
                              "-5.308229091527655e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 8   "
                              "-5.287106143354928e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8       2 9    "
                              "4.150843319352120e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 8      2 10   "
                              "-6.239469776425440e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 7       2 8       2 8   "
                              "-1.235972933403907e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 1   "
                              "-2.247776729083673e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 2    "
                              "6.536053782709580e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 3   "
                              "-2.394602011086554e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 4   "
                              "-1.009548657688206e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 5   "
                              "-4.610877925549496e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 6   "
                              "-1.035250899127595e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 7    "
                              "6.825000713784974e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 8    "
                              "6.793206654107093e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8       2 9   "
                              "-5.339205492488950e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 8      2 10    "
                              "7.986694351292455e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 8       2 6    "
                              "1.448191594496957e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 8       2 8   "
                              "-9.457805337346691e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 8      2 10   "
                              "-1.123688547349832e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8       2 2   "
                              "-1.998170202978424e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8       2 3    "
                              "7.486746703671627e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8       2 4    "
                              "3.108631840358286e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8       2 5    "
                              "1.533577798371569e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8       2 6    "
                              "3.244412860870672e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8       2 7   "
                              "-2.120590029749823e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8       2 8   "
                              "-2.118026380640986e+03 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8       2 9    "
                              "1.655194374023424e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 8      2 10   "
                              "-2.518089954832010e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 9       2 6   "
                              "-1.070248767816951e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 9       2 7    "
                              "4.253366829644002e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 9       2 8   "
                              "-4.552211551146173e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 9       2 9   "
                              "-4.633511357968799e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0       2 9      2 10   "
                              "-1.979024442684463e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 9       2 4    "
                              "1.016320947901380e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 9       2 5   "
                              "-3.039926879078884e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 9       2 6    "
                              "3.403319184508233e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 9       2 7   "
                              "-1.363257048405140e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 9       2 8    "
                              "1.452317475658280e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 9       2 9    "
                              "1.478825902886419e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 2       2 9      2 10    "
                              "6.329024260738756e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 9       2 3   "
                              "-1.395412567391564e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 9       2 4    "
                              "4.282577880577527e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 9       2 5   "
                              "-1.291638838655701e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 9       2 6    "
                              "1.436182314881393e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 9       2 7   "
                              "-5.825026979372014e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 9       2 8    "
                              "6.161858483127648e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 9       2 9    "
                              "6.278283290908753e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4       2 9      2 10    "
                              "2.696606804509701e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9       2 2   "
                              "-3.303255145773645e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9       2 3   "
                              "-9.637227239216637e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9       2 4    "
                              "2.942090964056198e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9       2 5   "
                              "-8.876029555580890e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9       2 6    "
                              "9.902405437257963e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9       2 7   "
                              "-3.924516205341249e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9       2 8    "
                              "4.150843319352119e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9       2 9    "
                              "4.273668554975398e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6       2 9      2 10    "
                              "1.792916316377662e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9       2 2    "
                              "4.261223997535471e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9       2 3    "
                              "1.243370880031008e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9       2 4   "
                              "-3.794591879474157e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9       2 5    "
                              "1.144451917307040e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9       2 6   "
                              "-1.277659427864181e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9       2 7    "
                              "5.046870394011087e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9       2 8   "
                              "-5.339205492488950e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9       2 9   "
                              "-5.504217097371161e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 8       2 9      2 10   "
                              "-2.302043226519036e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 9       2 9       2 9    "
                              "7.590235680273143e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9       2 2   "
                              "-1.309575084929303e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9       2 3   "
                              "-3.798183267541390e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9       2 4    "
                              "1.165758245058045e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9       2 5   "
                              "-3.505700739697208e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9       2 6    "
                              "3.915235751316376e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9       2 7   "
                              "-1.560969804035066e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9       2 8    "
                              "1.655194374023424e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9       2 9    "
                              "1.695819421675624e+03 \n";
    integralFileTwoBodyFAD << "       1 8      1 10       2 9      2 10    "
                              "7.180694203014171e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 0      2 10       2 6    "
                              "1.622634053486751e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0      2 10       2 8    "
                              "7.001489898205731e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0      2 10       2 9   "
                              "-1.979024442684464e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 0      2 10      2 10    "
                              "7.785748816323361e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2      2 10       2 4   "
                              "-1.255070378117838e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2      2 10       2 6   "
                              "-5.173506711398749e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2      2 10       2 7    "
                              "2.480996473328200e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2      2 10       2 8   "
                              "-2.238895790085912e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 2      2 10       2 9    "
                              "6.329024260738755e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 2      2 10      2 10   "
                              "-2.484293963254502e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 4      2 10       2 4   "
                              "-5.302150443587363e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4      2 10       2 5    "
                              "1.971049907245196e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 4      2 10       2 6   "
                              "-2.190682304157160e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4      2 10       2 7    "
                              "1.048492818872039e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4      2 10       2 8   "
                              "-9.591808561875582e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4      2 10       2 9    "
                              "2.696606804509701e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 4      2 10      2 10   "
                              "-1.055920326425384e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10       2 2    "
                              "3.372852130492381e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10       2 3   "
                              "-5.348072594774358e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10       2 4   "
                              "-3.567575078941885e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10       2 5    "
                              "1.391957884620223e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10       2 6   "
                              "-1.459885354939898e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10       2 7    "
                              "7.114493677598679e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10       2 8   "
                              "-6.239469776425442e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10       2 9    "
                              "1.792916316377662e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 6      2 10      2 10   "
                              "-6.997127383543939e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 7      2 10      2 10   "
                              "-1.653820579459890e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10       2 2   "
                              "-4.345238508775923e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10       2 3    "
                              "6.849785581267599e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10       2 4    "
                              "4.589079566073315e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10       2 5   "
                              "-1.801466221487417e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10       2 6    "
                              "1.875446118161311e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10       2 7   "
                              "-9.161076408066576e+01 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10       2 8    "
                              "7.986694351292452e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10       2 9   "
                              "-2.302043226519035e+02 \n";
    integralFileTwoBodyFAD << "       1 8       1 8      2 10      2 10    "
                              "8.980978031369179e+03 \n";
    integralFileTwoBodyFAD << "       1 8       1 9      2 10       2 8   "
                              "-1.123688547349832e+00 \n";
    integralFileTwoBodyFAD << "       1 8       1 9      2 10      2 10   "
                              "-1.254030655428824e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10       2 2    "
                              "1.342331337795869e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10       2 3   "
                              "-2.166530796206870e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10       2 4   "
                              "-1.425656523948271e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10       2 5    "
                              "5.454186991059916e+00 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10       2 6   "
                              "-5.854545589589042e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10       2 7    "
                              "2.832041378888439e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10       2 8   "
                              "-2.518089954832009e+02 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10       2 9    "
                              "7.180694203014171e+01 \n";
    integralFileTwoBodyFAD << "       1 8      1 10      2 10      2 10   "
                              "-2.808911020908904e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 0       2 1   "
                              "-1.730807465470137e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 0       2 1   "
                              "-6.555720287020161e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 0       2 1    "
                              "1.110179599027562e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 0       2 0    "
                              "1.799256995632455e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 0       2 1   "
                              "-3.357053977900662e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 0       2 2    "
                              "2.581993335835812e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 0       2 3    "
                              "1.961594360186819e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 0       2 1   "
                              "-1.609511939947150e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 0       2 0   "
                              "-2.319083748297069e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 0       2 1    "
                              "3.754632790792959e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 0       2 2   "
                              "-2.845826142019610e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 0       2 3   "
                              "-2.075563127712271e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 0       2 4   "
                              "-3.530267152205838e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 0       2 5    "
                              "1.724638751716351e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 0       2 7    "
                              "1.050604873883750e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 1       2 0   "
                              "-1.730807465470137e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 1       2 2    "
                              "2.451597067677715e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 1       2 0   "
                              "-6.555720287020161e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 1       2 2    "
                              "9.280444835313173e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 1       2 0    "
                              "1.110179599027562e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 1       2 2   "
                              "-1.572565038600371e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 1       2 3   "
                              "-1.588406858635854e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 1       2 0   "
                              "-3.357053977900662e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 1       2 1   "
                              "-2.137894636141125e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 1       2 2    "
                              "4.750875450135381e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 1       2 3    "
                              "4.395271864661643e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 1       2 4    "
                              "6.781248699377388e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 1       2 6    "
                              "1.435419026291401e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 1       2 0   "
                              "-1.609511939947150e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 1       2 2    "
                              "2.277692492269796e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 0    "
                              "3.754632790792959e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 1    "
                              "1.488789585354112e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 2   "
                              "-5.312282683849566e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 3   "
                              "-4.845690368846583e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 4   "
                              "-7.333908808192832e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 5    "
                              "9.708997917886078e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 6   "
                              "-1.610014872575666e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 7    "
                              "5.596731551590230e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 8   "
                              "-7.265476876791159e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1       2 9    "
                              "2.294151860610259e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 1      2 10    "
                              "3.208620506374662e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 2       2 1    "
                              "2.451597067677715e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 2       2 3   "
                              "-2.964881898903597e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 2       2 1    "
                              "9.280444835313173e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 2       2 3   "
                              "-1.121720773849247e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 2       2 4   "
                              "-1.340400183023981e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 2       2 1   "
                              "-1.572565038600372e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 2       2 3    "
                              "1.901626899321122e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 2       2 4    "
                              "2.525909890768940e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 0    "
                              "2.581993335835809e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 1    "
                              "4.750875450135381e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 2   "
                              "-6.583433655063617e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 3   "
                              "-5.740122322743413e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 4   "
                              "-7.050524343971145e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 5   "
                              "-1.756238495991596e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 7   "
                              "-2.448791263423931e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 8   "
                              "-1.861023223998612e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2       2 9   "
                              "-1.220413371137132e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 2      2 10    "
                              "1.251596804748747e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 2       2 1    "
                              "2.277692492269796e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 2       2 3   "
                              "-2.751898147320722e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 0   "
                              "-2.845826142019611e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 1   "
                              "-5.312282683849566e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 2    "
                              "5.914165703814945e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 3    "
                              "6.416940122247222e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 4    "
                              "7.782798621362857e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 5    "
                              "2.246132337274226e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 7    "
                              "2.740942830514093e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 8    "
                              "2.091920288798988e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2       2 9    "
                              "1.367956910899529e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 2      2 10   "
                              "-1.399146141323150e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 3       2 2   "
                              "-2.964881898903597e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 3       2 4   "
                              "-3.593929799520224e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 3       2 2   "
                              "-1.121720773849247e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 3       2 4   "
                              "-1.358501369865259e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 3       2 6    "
                              "1.000129368049821e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 3       2 1   "
                              "-1.588406858635854e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 3       2 2    "
                              "1.901626899321122e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 3       2 4    "
                              "2.307273346932322e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 3       2 5   "
                              "-1.183947176249873e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 3       2 6   "
                              "-1.820141439304317e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 3       2 4    "
                              "1.167903155640163e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 0    "
                              "1.961594360186821e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 1    "
                              "4.395271864661643e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 2   "
                              "-5.740122322743413e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 3    "
                              "1.066770086804321e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 4   "
                              "-6.952975158208280e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 5    "
                              "3.014409335018223e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 6    "
                              "5.147142479302496e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 7   "
                              "-1.043490069029153e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 8    "
                              "6.990737527840236e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3       2 9   "
                              "-3.536361341968031e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 3      2 10   "
                              "-2.024942526576488e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 3       2 2   "
                              "-2.751898147320722e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 3       2 4   "
                              "-3.333024663302270e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 0   "
                              "-2.075563127712262e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 1   "
                              "-4.845690368846592e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 2    "
                              "6.416940122247222e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 3   "
                              "-1.399495461472266e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 4    "
                              "7.770349593196544e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 5   "
                              "-3.276772253961199e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 6   "
                              "-5.679651405678728e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 7    "
                              "1.175838001351463e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 8   "
                              "-7.763585537812764e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3       2 9    "
                              "3.978505941123953e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 3      2 10    "
                              "2.236171432103814e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 4       2 3   "
                              "-3.593929799520224e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 4       2 5    "
                              "3.114070807786515e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 4       2 6    "
                              "1.247838238673187e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 4       2 2   "
                              "-1.340400183023981e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 4       2 3   "
                              "-1.358501369865259e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 4       2 4   "
                              "-2.328385456519772e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 4       2 5    "
                              "1.178916342433937e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 4       2 6    "
                              "4.625288618906221e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 4       2 2    "
                              "2.525909890768941e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 4       2 3    "
                              "2.307273346932322e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 4       2 4    "
                              "4.516657380965685e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 4       2 5   "
                              "-1.990417643706984e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 4       2 6   "
                              "-8.241155019424030e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 4       2 3    "
                              "1.167903155640163e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 4       2 5   "
                              "-1.010903912981626e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 1    "
                              "6.781248699377415e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 2   "
                              "-7.050524343971145e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 3   "
                              "-6.952975158208280e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 4   "
                              "-1.234468756169917e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 5    "
                              "6.019892800671735e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 6    "
                              "2.391210127227887e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 7    "
                              "2.846648138015005e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 8    "
                              "2.898027478814571e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4       2 9    "
                              "1.086284764807786e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 4      2 10   "
                              "-1.330019087231057e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 4       2 3   "
                              "-3.333024663302269e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 4       2 5    "
                              "2.886762933707562e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 4       2 6    "
                              "1.142896303026985e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 0   "
                              "-3.530267152205723e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 1   "
                              "-7.333908808192810e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 2    "
                              "7.782798621362859e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 3    "
                              "7.770349593196545e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 4    "
                              "1.357782188806382e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 5   "
                              "-6.730220282562677e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 6   "
                              "-2.655208117052649e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 7   "
                              "-3.207851805645732e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 8   "
                              "-3.244249965332596e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4       2 9   "
                              "-1.217963354263751e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 4      2 10    "
                              "1.482540158258932e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 5       2 4    "
                              "3.114070807786515e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 5       2 5    "
                              "2.799026088724463e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 5       2 6   "
                              "-3.819670485141559e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 5       2 4    "
                              "1.178916342433937e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 5       2 5    "
                              "1.075395185438978e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 5       2 6   "
                              "-1.439389276842473e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 5       2 7   "
                              "-2.768707399622547e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 5       2 3   "
                              "-1.183947176249874e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 5       2 4   "
                              "-1.990417643706984e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 5       2 5   "
                              "-1.745176364199179e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 5       2 6    "
                              "2.460688425725434e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 5       2 7    "
                              "4.460267096941017e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 5       2 9    "
                              "1.043158182951092e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 5       2 4   "
                              "-1.010903912981626e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 5       2 6    "
                              "1.239387730666073e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5       2 2   "
                              "-1.756238495991599e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5       2 3    "
                              "3.014409335018219e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5       2 4    "
                              "6.019892800671734e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5       2 5    "
                              "5.441386129562740e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5       2 6   "
                              "-7.372631609575187e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5       2 7   "
                              "-1.391514447957063e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5       2 8    "
                              "1.439431409543737e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5       2 9   "
                              "-3.265023418343503e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 5      2 10    "
                              "5.074491052839459e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 5       2 4    "
                              "2.886762933707562e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 5       2 5    "
                              "2.615391789411825e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 5       2 6   "
                              "-3.532842787495773e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 0    "
                              "1.724638751716437e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 1    "
                              "9.708997917886078e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 2    "
                              "2.246132337274226e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 3   "
                              "-3.276772253961204e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 4   "
                              "-6.730220282562677e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 5   "
                              "-6.112287151132913e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 6    "
                              "8.230527582700886e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 7    "
                              "1.561681332485156e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 8   "
                              "-1.551547624889774e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5       2 9    "
                              "3.667826103158186e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 5      2 10   "
                              "-5.733853449056491e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 6       2 4    "
                              "1.247838238673187e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 6       2 5   "
                              "-3.819670485141559e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 6       2 6   "
                              "-5.198931969245061e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 6       2 7    "
                              "1.647355097856701e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 6       2 8    "
                              "1.591502419242153e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 6       2 3    "
                              "1.000129368049822e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 6       2 4    "
                              "4.625288618906221e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 6       2 5   "
                              "-1.439389276842473e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 6       2 6   "
                              "-1.943779318542299e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 6       2 7    "
                              "6.239953534808498e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 6       2 8    "
                              "5.897004744747068e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 6      2 10   "
                              "-1.060383707821511e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 6       2 3   "
                              "-1.820141439304317e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 6       2 4   "
                              "-8.241155019424033e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 6       2 5    "
                              "2.460688425725434e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 6       2 6    "
                              "3.393520885342536e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 6       2 7   "
                              "-1.049786435114060e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 6       2 8   "
                              "-1.041815763577073e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 6       2 9   "
                              "-1.177377020620833e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 6      2 10    "
                              "1.866941760448747e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 6       2 5    "
                              "1.239387730666073e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 6       2 6    "
                              "1.685541059223235e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6       2 1    "
                              "1.435419026291401e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6       2 3    "
                              "5.147142479302490e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6       2 4    "
                              "2.391210127227887e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6       2 5   "
                              "-7.372631609575186e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6       2 6   "
                              "-1.000866563941849e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6       2 7    "
                              "3.180067623445261e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6       2 8    "
                              "3.029150003666393e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6       2 9    "
                              "3.647176483226015e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 6      2 10   "
                              "-5.464381642888395e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 6       2 4    "
                              "1.142896303026985e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 6       2 5   "
                              "-3.532842787495773e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 6       2 6   "
                              "-4.789834520304604e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 6       2 7    "
                              "1.525420447878576e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 6       2 8    "
                              "1.448191594496759e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6       2 1   "
                              "-1.610014872575712e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6       2 3   "
                              "-5.679651405678728e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6       2 4   "
                              "-2.655208117052649e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6       2 5    "
                              "8.230527582700886e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6       2 6    "
                              "1.114551578744252e+04 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6       2 7   "
                              "-3.555416616107038e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6       2 8   "
                              "-3.361825765222610e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6       2 9   "
                              "-4.095087705186355e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 6      2 10    "
                              "6.075992096087928e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 7       2 6    "
                              "1.647355097856701e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 7       2 7    "
                              "6.456504337567213e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 7       2 8   "
                              "-1.029351659267908e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 7       2 5   "
                              "-2.768707399622547e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 7       2 6    "
                              "6.239953534808498e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 7       2 7    "
                              "2.462404879840944e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 7       2 8   "
                              "-3.864053571595490e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 7       2 9   "
                              "-2.882133528804211e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 7       2 5    "
                              "4.460267096941013e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 7       2 6   "
                              "-1.049786435114060e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 7       2 7   "
                              "-4.070465056235976e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 7       2 8    "
                              "6.645736996135041e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 7       2 9    "
                              "4.774751388295488e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 7       2 7   "
                              "-2.096101407636390e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7       2 2   "
                              "-2.448791263423925e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7       2 3   "
                              "-1.043490069029153e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7       2 4    "
                              "2.846648138015011e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7       2 5   "
                              "-1.391514447957064e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7       2 6    "
                              "3.180067623445261e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7       2 7    "
                              "1.250162959060689e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7       2 8   "
                              "-1.978191421576763e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7       2 9   "
                              "-1.455209600846863e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 7      2 10    "
                              "2.640636332093961e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 7       2 6    "
                              "1.525420447878576e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 7       2 7    "
                              "6.003160326643294e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 0    "
                              "1.050604873883779e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 1    "
                              "5.596731551590230e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 2    "
                              "2.740942830514122e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 3    "
                              "1.175838001351461e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 4   "
                              "-3.207851805645732e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 5    "
                              "1.561681332485156e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 6   "
                              "-3.555416616107038e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 7   "
                              "-1.400819391069174e+04 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 8    "
                              "2.205223921682365e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7       2 9    "
                              "1.626725811013738e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 7      2 10   "
                              "-2.951224092558840e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 8       2 6    "
                              "1.591502419242153e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 8       2 7   "
                              "-1.029351659267908e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 8       2 8   "
                              "-1.030555941314443e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 8      2 10   "
                              "-1.252907494784057e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 8       2 6    "
                              "5.897004744747069e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 8       2 7   "
                              "-3.864053571595490e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 8       2 8   "
                              "-3.845511104588639e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 8       2 9    "
                              "3.023552700596141e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 8      2 10   "
                              "-4.610521400014492e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 8       2 6   "
                              "-1.041815763577073e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 8       2 7    "
                              "6.645736996135041e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 8       2 8    "
                              "6.718945757846795e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 8       2 9   "
                              "-5.144945807610154e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 8      2 10    "
                              "8.157453662324455e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 8       2 8    "
                              "3.330593692195275e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8       2 2   "
                              "-1.861023223998613e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8       2 3    "
                              "6.990737527840220e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8       2 4    "
                              "2.898027478814566e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8       2 5    "
                              "1.439431409543741e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8       2 6    "
                              "3.029150003666393e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8       2 7   "
                              "-1.978191421576763e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8       2 8   "
                              "-1.976566944719782e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8       2 9    "
                              "1.543651779489834e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 8      2 10   "
                              "-2.352073556074061e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 8       2 6    "
                              "1.448191594496758e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 8       2 8   "
                              "-9.457805337346395e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 8      2 10   "
                              "-1.123688547349939e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 1   "
                              "-7.265476876791230e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 2    "
                              "2.091920288798976e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 3   "
                              "-7.763585537812742e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 4   "
                              "-3.244249965332597e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 5   "
                              "-1.551547624889774e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 6   "
                              "-3.361825765222611e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 7    "
                              "2.205223921682365e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 8    "
                              "2.199302889899705e+04 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8       2 9   "
                              "-1.722930599907308e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 8      2 10    "
                              "2.602987225765866e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 1       2 9       2 9    "
                              "8.178121829989717e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 9       2 7   "
                              "-2.882133528804211e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 9       2 8    "
                              "3.023552700596140e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 9       2 9    "
                              "3.115214987854688e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3       2 9      2 10    "
                              "1.313011779560913e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 9       2 5    "
                              "1.043158182951092e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 9       2 6   "
                              "-1.177377020620833e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 9       2 7    "
                              "4.774751388295488e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 9       2 8   "
                              "-5.144945807610155e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 9       2 9   "
                              "-5.157927070990905e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5       2 9      2 10   "
                              "-2.266214214577090e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 6       2 9       2 9   "
                              "-2.651444217584577e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9       2 2   "
                              "-1.220413371137132e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9       2 3   "
                              "-3.536361341968040e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9       2 4    "
                              "1.086284764807789e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9       2 5   "
                              "-3.265023418343503e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9       2 6    "
                              "3.647176483226014e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9       2 7   "
                              "-1.455209600846863e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9       2 8    "
                              "1.543651779489834e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9       2 9    "
                              "1.580454111766436e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 7       2 9      2 10    "
                              "6.700652879212963e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 8       2 9       2 9    "
                              "7.590235680273349e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 1    "
                              "2.294151860610245e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 2    "
                              "1.367956910899522e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 3    "
                              "3.978505941123942e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 4   "
                              "-1.217963354263753e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 5    "
                              "3.667826103158184e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 6   "
                              "-4.095087705186355e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 7    "
                              "1.626725811013738e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 8   "
                              "-1.722930599907307e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9       2 9   "
                              "-1.769898187824457e+04 \n";
    integralFileTwoBodyFAD << "       1 9       1 9       2 9      2 10   "
                              "-7.455773271252704e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 1      2 10       2 8   "
                              "-1.252907494784057e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 1      2 10      2 10   "
                              "-1.375296255861436e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 3      2 10       2 6   "
                              "-1.060383707821511e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3      2 10       2 8   "
                              "-4.610521400014492e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3      2 10       2 9    "
                              "1.313011779560913e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 3      2 10      2 10   "
                              "-5.111371542171212e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 5      2 10       2 6    "
                              "1.866941760448747e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5      2 10       2 8    "
                              "8.157453662324455e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5      2 10       2 9   "
                              "-2.266214214577090e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 5      2 10      2 10    "
                              "8.964759125910372e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 6      2 10      2 10    "
                              "4.422719749063024e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10       2 2    "
                              "1.251596804748747e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10       2 3   "
                              "-2.024942526576488e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10       2 4   "
                              "-1.330019087231056e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10       2 5    "
                              "5.074491052839461e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10       2 6   "
                              "-5.464381642888395e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10       2 7    "
                              "2.640636332093961e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10       2 8   "
                              "-2.352073556074062e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10       2 9    "
                              "6.700652879212963e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 7      2 10      2 10   "
                              "-2.621992222462075e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 8      2 10       2 8   "
                              "-1.123688547349939e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 8      2 10      2 10   "
                              "-1.254030655428823e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 1    "
                              "3.208620506374705e+00 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 2   "
                              "-1.399146141323150e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 3    "
                              "2.236171432103808e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 4    "
                              "1.482540158258932e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 5   "
                              "-5.733853449056480e+01 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 6    "
                              "6.075992096087925e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 7   "
                              "-2.951224092558841e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 8    "
                              "2.602987225765865e+03 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10       2 9   "
                              "-7.455773271252704e+02 \n";
    integralFileTwoBodyFAD << "       1 9       1 9      2 10      2 10    "
                              "2.913012974614897e+04 \n";
    integralFileTwoBodyFAD << "       1 9      1 10      2 10      2 10   "
                              "-1.207042296544339e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 0       2 1   "
                              "-3.125876669004099e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 0       2 1   "
                              "-5.718389350457165e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 0       2 1    "
                              "1.172883015425946e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 0       2 2   "
                              "-9.142470620259783e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 0       2 1    "
                              "1.503299029713919e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 0       2 0    "
                              "1.992058670279703e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 0       2 1   "
                              "-3.600904551175626e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 0       2 2    "
                              "2.761058304293006e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 0       2 3    "
                              "2.082722879848016e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 0       2 0   "
                              "-2.313025647662081e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 0       2 1    "
                              "3.743835247339773e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 0       2 2   "
                              "-2.837558014729725e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 0       2 3   "
                              "-2.069336045400967e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 0       2 4   "
                              "-3.519782892750927e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 0       2 5    "
                              "1.719679356825350e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 0       2 7    "
                              "1.047587345595131e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 1       2 2   "
                              "-1.292872672988605e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 1       2 0   "
                              "-3.125876669004100e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 1       2 2    "
                              "4.426763503373087e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 1       2 0   "
                              "-5.718389350457165e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 1       2 2    "
                              "8.094570959130648e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 1       2 0    "
                              "1.172883015425946e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 1       2 1    "
                              "1.008383095100763e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 1       2 2   "
                              "-1.660150637030747e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 1       2 3   "
                              "-1.555830872314816e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 1       2 4   "
                              "-2.434022040650373e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 1       2 0    "
                              "1.503299029713919e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 1       2 2   "
                              "-2.127438457360219e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 1       2 0   "
                              "-3.600904551175626e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 1       2 1   "
                              "-2.110613471131293e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 1       2 2    "
                              "5.095763963044434e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 1       2 3    "
                              "4.700424631374045e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 1       2 4    "
                              "7.228580193993204e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 1       2 6    "
                              "1.539879421310991e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 0    "
                              "3.743835247339773e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 1    "
                              "1.482702850793540e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 2   "
                              "-5.297002808300513e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 3   "
                              "-4.831612794067354e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 4   "
                              "-7.312270423654239e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 5    "
                              "9.680956501355853e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 6   "
                              "-1.605405266754068e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 7    "
                              "5.580740231146485e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 8   "
                              "-7.244527851895711e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1       2 9    "
                              "2.287591315331744e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 1      2 10    "
                              "3.199337758092716e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 2       2 1   "
                              "-1.292872672988605e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 2       2 3    "
                              "1.563931832604049e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 2       2 1    "
                              "4.426763503373087e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 2       2 3   "
                              "-5.352581309875099e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 2       2 1    "
                              "8.094570959130648e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 2       2 3   "
                              "-9.783918468214852e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 2       2 4   "
                              "-1.120380157150106e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 2       2 0   "
                              "-9.142470620259790e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 2       2 1   "
                              "-1.660150637030747e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 2       2 2    "
                              "2.719105181960352e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 2       2 3    "
                              "2.006180858615087e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 2       2 4    "
                              "2.492695068473069e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 2       2 1   "
                              "-2.127438457360219e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 2       2 3    "
                              "2.570381918532745e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 0    "
                              "2.761058304293001e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 1    "
                              "5.095763963044434e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 2   "
                              "-6.769114163844168e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 3   "
                              "-6.156582684296372e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 4   "
                              "-7.542149448359859e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 5   "
                              "-1.934359721557095e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 7   "
                              "-2.626784858720399e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 8   "
                              "-1.998170202978413e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2       2 9   "
                              "-1.309575084929275e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 2      2 10    "
                              "1.342331337795865e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 0   "
                              "-2.837558014729726e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 1   "
                              "-5.297002808300513e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 2    "
                              "5.894236859758028e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 3    "
                              "6.398479435195523e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 4    "
                              "7.760205813246472e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 5    "
                              "2.240283009396094e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 7    "
                              "2.733066952364615e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 8    "
                              "2.085925927748521e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2       2 9    "
                              "1.364029692171971e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 2      2 10   "
                              "-1.395121494182030e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 3       2 2    "
                              "1.563931832604049e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 3       2 4    "
                              "1.896344383368161e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 3       2 2   "
                              "-5.352581309875100e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 3       2 4   "
                              "-6.486161216250300e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 3       2 2   "
                              "-9.783918468214853e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 3       2 4   "
                              "-1.184225594763195e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 1   "
                              "-1.555830872314816e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 2    "
                              "2.006180858615087e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 3   "
                              "-3.131303957664490e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 4    "
                              "2.430724079773960e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 5   "
                              "-1.080761725485803e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 6   "
                              "-1.818513046857072e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 7    "
                              "3.618890254776677e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 8   "
                              "-2.459244785532984e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 3       2 9    "
                              "1.228450458146236e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 3       2 2    "
                              "2.570381918532745e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 3       2 4    "
                              "3.113584990697563e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 0    "
                              "2.082722879848023e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 1    "
                              "4.700424631374045e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 2   "
                              "-6.156582684296372e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 3    "
                              "1.185843945031393e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 4   "
                              "-7.456977486817236e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 5    "
                              "3.214129084211634e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 6    "
                              "5.506918924840444e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 7   "
                              "-1.121168330839431e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 8    "
                              "7.486746703671661e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3       2 9   "
                              "-3.798183267541409e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 3      2 10   "
                              "-2.166530796206880e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 0   "
                              "-2.069336045400967e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 1   "
                              "-4.831612794067350e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 2    "
                              "6.398479435195523e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 3   "
                              "-1.395889857703450e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 4    "
                              "7.747990076349642e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 5   "
                              "-3.267153571552591e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 6   "
                              "-5.663148293082413e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 7    "
                              "1.172471457256250e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 8   "
                              "-7.741163077347272e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3       2 9    "
                              "3.967108432514041e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 3      2 10    "
                              "2.229680573063598e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 4       2 3    "
                              "1.896344383368161e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 4       2 5   "
                              "-1.642475989091995e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 4       2 3   "
                              "-6.486161216250300e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 4       2 4   "
                              "-1.147919633989249e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 4       2 5    "
                              "5.623506490270193e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 4       2 6    "
                              "2.236000088444868e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 4       2 2   "
                              "-1.120380157150106e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 4       2 3   "
                              "-1.184225594763195e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 4       2 4   "
                              "-1.921375436823263e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 4       2 5    "
                              "1.030261558620388e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 4       2 6    "
                              "3.960782823022988e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 1   "
                              "-2.434022040650373e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 2    "
                              "2.492695068473068e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 3    "
                              "2.430724079773960e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 4    "
                              "4.378540156496537e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 5   "
                              "-2.103620981509875e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 6   "
                              "-8.407319943652762e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 7   "
                              "-9.883299497889421e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 8   "
                              "-1.012354241756134e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4       2 9   "
                              "-3.787129376070280e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 4      2 10    "
                              "4.660740009072184e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 4       2 3    "
                              "3.113584990697563e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 4       2 5   "
                              "-2.695151589885695e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 4       2 6   "
                              "-1.072252078131462e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 1    "
                              "7.228580193993204e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 2   "
                              "-7.542149448359859e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 3   "
                              "-7.456977486817236e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 4   "
                              "-1.319562657339694e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 5    "
                              "6.456887602996155e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 6    "
                              "2.561210028161615e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 7    "
                              "3.057786407931134e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 8    "
                              "3.108631840358292e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4       2 9    "
                              "1.165758245058052e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 4      2 10   "
                              "-1.425656523948273e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 0   "
                              "-3.519782892750928e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 1   "
                              "-7.312270423654239e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 2    "
                              "7.760205813246470e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 3    "
                              "7.747990076349642e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 4    "
                              "1.353830425661771e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 5   "
                              "-6.710858557077943e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 6   "
                              "-2.647531339312153e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 7   "
                              "-3.198675616300079e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 8   "
                              "-3.234933504581290e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4       2 9   "
                              "-1.214466982556261e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 4      2 10    "
                              "1.478266260459662e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 5       2 4   "
                              "-1.642475989091995e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 5       2 5   "
                              "-1.468996738069969e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 5       2 6    "
                              "2.017684486939139e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 5       2 4    "
                              "5.623506490270193e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 5       2 5    "
                              "5.082349703374113e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 5       2 6   "
                              "-6.885951321223107e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 5       2 7   "
                              "-1.309457270819691e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 5       2 4    "
                              "1.030261558620388e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 5       2 5    "
                              "9.531551167504970e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 5       2 6   "
                              "-1.252039584713829e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 5       2 7   "
                              "-2.464810671566426e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 5       2 3   "
                              "-1.080761725485801e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 5       2 4   "
                              "-2.103620981509875e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 5       2 5   "
                              "-1.893254901229216e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 5       2 6    "
                              "2.579789406880491e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 5       2 7    "
                              "4.843737105543599e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 5       2 8   "
                              "-5.178464362429572e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 5       2 9    "
                              "1.135705037198188e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 5      2 10   "
                              "-1.757020148799411e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 5       2 4   "
                              "-2.695151589885695e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 5       2 5   "
                              "-2.433439865912232e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 5       2 6    "
                              "3.301975296170601e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5       2 2   "
                              "-1.934359721557095e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5       2 3    "
                              "3.214129084211637e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5       2 4    "
                              "6.456887602996155e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5       2 5    "
                              "5.842113970981339e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5       2 6   "
                              "-7.905409653530523e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5       2 7   "
                              "-1.493836520386779e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5       2 8    "
                              "1.533577798371563e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5       2 9   "
                              "-3.505700739697201e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 5      2 10    "
                              "5.454186991059953e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 0    "
                              "1.719679356825464e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 1    "
                              "9.680956501355498e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 2    "
                              "2.240283009396072e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 3   "
                              "-3.267153571552591e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 4   "
                              "-6.710858557077943e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 5   "
                              "-6.094763108442174e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 6    "
                              "8.206824823600957e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 7    "
                              "1.557201326843723e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 8   "
                              "-1.546951545941612e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5       2 9    "
                              "3.657307585543252e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 5      2 10   "
                              "-5.717489191718133e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 6       2 5    "
                              "2.017684486939138e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 6       2 6    "
                              "2.753394733772482e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 6       2 4    "
                              "2.236000088444868e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 6       2 5   "
                              "-6.885951321223107e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 6       2 6   "
                              "-9.345209158285757e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 6       2 7    "
                              "2.975652051047463e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 6       2 8    "
                              "2.851784702994153e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 6       2 4    "
                              "3.960782823022988e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 6       2 5   "
                              "-1.252039584713829e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 6       2 6   "
                              "-1.677232895230524e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 6       2 7    "
                              "5.462478120731864e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 6       2 8    "
                              "5.074230704757856e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 6       2 3   "
                              "-1.818513046857073e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 6       2 4   "
                              "-8.407319943652756e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 6       2 5    "
                              "2.579789406880491e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 6       2 6    "
                              "3.510151522098444e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 6       2 7   "
                              "-1.111132929897306e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 6       2 8   "
                              "-1.065161737061374e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 6       2 9   "
                              "-1.269762453613437e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 6      2 10    "
                              "1.918822020574653e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 6       2 4   "
                              "-1.072252078131462e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 6       2 5    "
                              "3.301975296170601e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 6       2 6    "
                              "4.485429484789082e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 6       2 7   "
                              "-1.423405578005336e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 6       2 8   "
                              "-1.357942837917659e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6       2 1    "
                              "1.539879421311005e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6       2 3    "
                              "5.506918924840444e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6       2 4    "
                              "2.561210028161615e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6       2 5   "
                              "-7.905409653530525e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6       2 6   "
                              "-1.072637165170560e+03 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6       2 7    "
                              "3.411000185590887e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6       2 8    "
                              "3.244412860870671e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6       2 9    "
                              "3.915235751316374e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 6      2 10   "
                              "-5.854545589589033e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6       2 1   "
                              "-1.605405266754078e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6       2 3   "
                              "-5.663148293082417e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6       2 4   "
                              "-2.647531339312154e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6       2 5    "
                              "8.206824823600957e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6       2 6    "
                              "1.111336034446191e+04 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6       2 7   "
                              "-3.545188467164468e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6       2 8   "
                              "-3.352099300934417e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6       2 9   "
                              "-4.083342826350892e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 6      2 10    "
                              "6.058443275953132e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 7       2 7   "
                              "-3.396837755305310e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 7       2 5   "
                              "-1.309457270819690e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 7       2 6    "
                              "2.975652051047463e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 7       2 7    "
                              "1.169202560602879e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 7       2 8   "
                              "-1.853176606030542e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 7       2 9   "
                              "-1.372848328346638e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 7       2 5   "
                              "-2.464810671566426e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 7       2 6    "
                              "5.462478120731864e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 7       2 7    "
                              "2.169086033921917e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 7       2 8   "
                              "-3.355787621948809e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 7       2 9   "
                              "-2.542766261503388e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 7       2 3    "
                              "3.618890254776677e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 7       2 4   "
                              "-9.883299497889420e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 7       2 5    "
                              "4.843737105543599e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 7       2 6   "
                              "-1.111132929897306e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 7       2 7   "
                              "-4.359417594991876e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 7       2 8    "
                              "6.929917642364798e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 7       2 9    "
                              "5.083215609747409e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 7      2 10   "
                              "-9.231586122277207e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 7       2 6   "
                              "-1.423405578005336e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 7       2 7   "
                              "-5.592939369227919e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7       2 2   "
                              "-2.626784858720399e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7       2 3   "
                              "-1.121168330839429e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7       2 4    "
                              "3.057786407931137e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7       2 5   "
                              "-1.493836520386779e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7       2 6    "
                              "3.411000185590887e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7       2 7    "
                              "1.341556434415410e+03 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7       2 8   "
                              "-2.120590029749823e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7       2 9   "
                              "-1.560969804035064e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 7      2 10    "
                              "2.832041378888436e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 0    "
                              "1.047587345595049e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 1    "
                              "5.580740231146556e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 2    "
                              "2.733066952364638e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 3    "
                              "1.172471457256243e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 4   "
                              "-3.198675616300074e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 5    "
                              "1.557201326843723e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 6   "
                              "-3.545188467164468e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 7   "
                              "-1.396796090211070e+04 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 8    "
                              "2.198865985710166e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7       2 9    "
                              "1.622047255417133e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 7      2 10   "
                              "-2.942733394899640e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 8       2 8    "
                              "5.463276728836940e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 8       2 6    "
                              "2.851784702994153e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 8       2 7   "
                              "-1.853176606030542e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 8       2 8   "
                              "-1.851225908332972e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 8       2 9    "
                              "1.446434211493660e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 8      2 10   "
                              "-2.239840669708080e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 8       2 6    "
                              "5.074230704757855e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 8       2 7   "
                              "-3.355787621948809e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 8       2 8   "
                              "-3.319416623364565e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 8       2 9    "
                              "2.636607645253483e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 8      2 10   "
                              "-3.974403634290994e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 8       2 3   "
                              "-2.459244785532985e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 8       2 4   "
                              "-1.012354241756135e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 8       2 5   "
                              "-5.178464362429572e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 8       2 6   "
                              "-1.065161737061373e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 8       2 7    "
                              "6.929917642364799e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 8       2 8    "
                              "6.936027991558964e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 8       2 9   "
                              "-5.401560295211443e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 8      2 10    "
                              "8.287079967044906e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 8       2 6   "
                              "-1.357942837917659e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 8       2 8    "
                              "8.858111689187419e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 8      2 10    "
                              "1.054592973375522e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8       2 2   "
                              "-1.998170202978441e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8       2 3    "
                              "7.486746703671646e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8       2 4    "
                              "3.108631840358288e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8       2 5    "
                              "1.533577798371562e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8       2 6    "
                              "3.244412860870671e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8       2 7   "
                              "-2.120590029749823e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8       2 8   "
                              "-2.118026380640986e+03 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8       2 9    "
                              "1.655194374023424e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 8      2 10   "
                              "-2.518089954832012e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 1   "
                              "-7.244527851895697e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 2    "
                              "2.085925927748520e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 3   "
                              "-7.741163077347292e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 4   "
                              "-3.234933504581293e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 5   "
                              "-1.546951545941612e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 6   "
                              "-3.352099300934418e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 7    "
                              "2.198865985710166e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 8    "
                              "2.192953203938431e+04 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8       2 9   "
                              "-1.717967914830158e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 8      2 10    "
                              "2.595435963734375e+03 \n";
    integralFileTwoBodyFAD << "      1 10       1 0       2 9       2 9   "
                              "-4.303707474426996e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 9       2 7   "
                              "-1.372848328346639e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 9       2 8    "
                              "1.446434211493660e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2       2 9       2 9    "
                              "1.480443351584373e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 9       2 7   "
                              "-2.542766261503385e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 9       2 8    "
                              "2.636607645253482e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 9       2 9    "
                              "2.743726361158308e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4       2 9      2 10    "
                              "1.140490862655851e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 9       2 3    "
                              "1.228450458146236e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 9       2 4   "
                              "-3.787129376070281e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 9       2 5    "
                              "1.135705037198188e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 9       2 6   "
                              "-1.269762453613438e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 9       2 7    "
                              "5.083215609747409e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 9       2 8   "
                              "-5.401560295211443e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 9       2 9   "
                              "-5.513634585383528e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 6       2 9      2 10   "
                              "-2.350614110853461e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 7       2 9       2 9   "
                              "-7.069141743283677e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9       2 2   "
                              "-1.309575084929275e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9       2 3   "
                              "-3.798183267541432e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9       2 4    "
                              "1.165758245058051e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9       2 5   "
                              "-3.505700739697205e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9       2 6    "
                              "3.915235751316376e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9       2 7   "
                              "-1.560969804035064e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9       2 8    "
                              "1.655194374023424e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9       2 9    "
                              "1.695819421675624e+03 \n";
    integralFileTwoBodyFAD << "      1 10       1 8       2 9      2 10    "
                              "7.180694203014171e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 1    "
                              "2.287591315331758e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 2    "
                              "1.364029692171966e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 3    "
                              "3.967108432514052e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 4   "
                              "-1.214466982556258e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 5    "
                              "3.657307585543252e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 6   "
                              "-4.083342826350891e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 7    "
                              "1.622047255417132e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 8   "
                              "-1.717967914830158e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9       2 9   "
                              "-1.764813087629521e+04 \n";
    integralFileTwoBodyFAD << "      1 10      1 10       2 9      2 10   "
                              "-7.434249539143888e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 0      2 10      2 10    "
                              "7.304753787658044e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2      2 10       2 8   "
                              "-2.239840669708080e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 2      2 10      2 10   "
                              "-2.467054026491808e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 4      2 10       2 8   "
                              "-3.974403634290994e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4      2 10       2 9    "
                              "1.140490862655851e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 4      2 10      2 10   "
                              "-4.409889931930826e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6      2 10       2 4    "
                              "4.660740009072189e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6      2 10       2 5   "
                              "-1.757020148799411e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6      2 10       2 6    "
                              "1.918822020574652e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6      2 10       2 7   "
                              "-9.231586122277209e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 6      2 10       2 8    "
                              "8.287079967044906e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6      2 10       2 9   "
                              "-2.350614110853462e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 6      2 10      2 10    "
                              "9.211453258183601e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 7      2 10       2 8    "
                              "1.054592973375523e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 7      2 10      2 10    "
                              "1.175238565629062e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10       2 2    "
                              "1.342331337795871e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10       2 3   "
                              "-2.166530796206874e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10       2 4   "
                              "-1.425656523948272e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10       2 5    "
                              "5.454186991059938e+00 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10       2 6   "
                              "-5.854545589589033e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10       2 7    "
                              "2.832041378888438e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10       2 8   "
                              "-2.518089954832012e+02 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10       2 9    "
                              "7.180694203014171e+01 \n";
    integralFileTwoBodyFAD << "      1 10       1 8      2 10      2 10   "
                              "-2.808911020908904e+03 \n";
    integralFileTwoBodyFAD << "      1 10       1 9      2 10      2 10   "
                              "-1.207042296543074e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 1    "
                              "3.199337758092726e+00 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 2   "
                              "-1.395121494182030e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 3    "
                              "2.229680573063597e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 4    "
                              "1.478266260459662e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 5   "
                              "-5.717489191718111e+01 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 6    "
                              "6.058443275953133e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 7   "
                              "-2.942733394899641e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 8    "
                              "2.595435963734375e+03 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10       2 9   "
                              "-7.434249539143888e+02 \n";
    integralFileTwoBodyFAD << "      1 10      1 10      2 10      2 10    "
                              "2.904591355839877e+04 \n";
    integralFileTwoBodyFAD.close();
    // Integral file for the fingerprint region
    integralFileTwoBodyFADFingerPrint.open(
        "integral_file_test_TwoBodyFAD_Fingerprint"
    );
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0   6.275787578691857e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1   1.494468129036738e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2   2.321175402307838e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3   2.117054739827139e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0   2.044645816296873e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1   1.881219476480249e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2  -1.699890107216312e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3  -7.224904760154720e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0   1.886164154398152e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1  -1.860509644906042e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2   3.136374492379533e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3  -9.518194826058090e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0   1.633708908816377e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1  -6.302829636380342e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2  -1.368240506245850e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3   4.406882178014263e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0   6.309851511299568e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1   1.086133852953992e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2  -1.962848829510437e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3  -1.755474877782550e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0   7.132139541049192e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1   1.896696099442117e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2  -1.339611560870598e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3   1.243314199153838e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0  -2.875907694514036e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1  -2.021872791566559e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2   3.171041668461655e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3   2.168363758811359e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0  -1.849541412586952e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1   1.430703659791857e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2   2.626426884815404e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3   4.461143405644639e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0   7.031334063379517e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1   3.709865155250298e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2  -1.968518847912678e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3   2.063143089445170e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0   2.783932559939969e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1   2.110577451071812e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2  -2.563220122298612e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3   2.751353191907417e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0  -2.447217122856030e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1  -1.503577338928864e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2   3.521815171401519e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3  -1.052950302632676e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0   1.447806666119930e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1   5.547322922061981e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2  -6.484270225252075e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3   4.945908029189794e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0   7.054220685870023e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1  -3.886751139175968e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2  -3.985744083411826e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3  -2.509680276402036e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0   7.773037522745194e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1   2.118194638827003e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2   1.264878859306626e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3  -8.001232525655082e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0  -1.015398049607994e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1   8.484578356665845e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2   3.536338932892386e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3   4.402165107051867e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0  -5.090003405028722e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1   6.856074710644969e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2   1.677581838092739e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3   4.968765472391655e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-0       5-0   7.241888750388754e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-0       5-1  -5.493876027286752e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-0       5-2  -1.007590155414743e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-0       5-3   1.137721725302256e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-1       5-0  -4.876215825927807e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-1       5-1   2.173330564227885e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-1       5-2   2.758789281615275e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-1       5-3   4.344163607330576e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-2       5-0  -3.119439632446382e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-2       5-1  -8.618427971407287e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-2       5-2   3.625571898722256e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-2       5-3  -2.718634397192842e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-3       5-0  -8.172262588935814e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-3       5-1   5.499731947931691e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-3       5-2  -1.397984865077282e-12 \n";
    integralFileTwoBodyFADFingerPrint
        << "5-3       5-3   5.090521322065309e+02 \n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-0       2-0    1.297343124233882e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-0       2-1    2.120292083830328e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-0       2-2    1.785012427390675e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-0       2-3   -7.288098960461179e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-0       2-0   -3.442142023426560e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-0       2-1   -9.932297685313982e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-0       2-2   -4.876555842404728e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-0       2-3   -1.150024365761365e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-0       2-0    5.603508807100630e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-0       2-1    3.011990983997775e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-0       2-2    7.865604329382982e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-0       2-3    1.168101584359989e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-0       2-0    5.803069401295355e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-0       2-1    8.783173312640950e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-0       2-2    8.136045835280147e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-0       2-3    1.000086050419885e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-1       2-0    2.120292083830264e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-1       2-1    3.808864440106895e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-1       2-2    3.082815193203547e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-1       2-3   -2.966939908014025e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-1       2-0   -9.932297685225164e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-1       2-1   -1.030366499686407e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-1       2-2   -1.814846013816568e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-1       2-3    8.422754729641764e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-1       2-0    3.011990983997884e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-1       2-1    1.667065778132682e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-1       2-2    4.422417099641768e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-1       2-3   -1.342531342686826e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-1       2-0    8.783173312640265e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-1       2-1    1.725072621922416e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-1       2-2    1.258453594400861e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-1       2-3   -1.386618534247728e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-2       2-0    1.785012427390675e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-2       2-1    3.082815193203654e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-2       2-2    6.140410507941815e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-2       2-3   -3.841030999845806e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-2       2-0   -4.876555842404728e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-2       2-1   -1.814846013816443e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-2       2-2   -1.705959768527694e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-2       2-3    2.294550993594324e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-2       2-0    7.865604329382982e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-2       2-1    4.422417099641768e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-2       2-2    2.737263975189050e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-2       2-3   -5.515976450455609e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-2       2-0    8.136045835280147e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-2       2-1    1.258453594400972e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-2       2-2    2.829543743995078e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-2       2-3   -1.547260227533053e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-3       2-0   -7.288098960371190e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-3       2-1   -2.966939908014025e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-3       2-2   -3.841030999845698e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       2-3       2-3    8.215221415257032e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-3       2-0   -1.150024365757805e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-3       2-1    8.422754729641765e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-3       2-2    2.294550993597877e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       2-3       2-3   -2.342082778975014e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-3       2-0    1.168101584358861e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-3       2-1   -1.342531342686826e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-3       2-2   -5.515976450455179e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       2-3       2-3    3.728313871997137e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-3       2-0    1.000086050421272e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-3       2-1   -1.386618534247728e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-3       2-2   -1.547260227533329e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       2-3       2-3    3.850284070723647e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-0       2-0   -3.442142023426559e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-0       2-1   -9.932297685296219e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-0       2-2   -4.876555842404728e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-0       2-3   -1.150024365764918e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-0       2-0    3.864817720631761e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-0       2-1    6.231093470653708e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-0       2-2    5.318974633135891e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-0       2-3   -2.125335204897823e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-0       2-0   -4.813650081948035e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-0       2-1   -2.036304517093291e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-0       2-2   -6.820273066847479e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-0       2-3   -1.612743304280462e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-0       2-0   -9.386389032235362e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-0       2-1   -5.127844970083305e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-0       2-2   -1.317634215388623e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-0       2-3   -1.948199092011351e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-1       2-0   -9.932297685225164e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-1       2-1   -1.030366499686407e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-1       2-2   -1.814846013802357e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-1       2-3    8.422754729641764e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-1       2-0    6.231093470653153e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-1       2-1    1.134864311043026e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-1       2-2    9.063039458064095e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-1       2-3   -8.843939648592936e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-1       2-0   -2.036304517093168e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-1       2-1   -1.441006732780512e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-1       2-2   -3.441438905755517e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-1       2-3    1.178141724896658e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-1       2-0   -5.127844970082805e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-1       2-1   -2.792591155995216e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-1       2-2   -7.525499252658753e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-1       2-3    2.249145748800627e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-2       2-0   -4.876555842404728e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-2       2-1   -1.814846013802232e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-2       2-2   -1.705959768527694e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-2       2-3    2.294550993594324e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-2       2-0    5.318974633135890e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-2       2-1    9.063039458065665e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-2       2-2    1.829990400655880e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-2       2-3   -1.129526381714806e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-2       2-0   -6.820273066847479e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-2       2-1   -3.441438905755517e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-2       2-2   -2.386061760883756e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-2       2-3    4.297719780578690e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-2       2-0   -1.317634215388623e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-2       2-1   -7.525499252658836e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-2       2-2   -4.585562739975720e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-2       2-3    9.384490959353176e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-3       2-0   -1.150024365761358e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-3       2-1    8.422754729641765e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-3       2-2    2.294550993597877e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       2-3       2-3   -2.342082778975015e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-3       2-0   -2.125335204928355e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-3       2-1   -8.843939648592935e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-3       2-2   -1.129526381714761e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       2-3       2-3    2.448900423435290e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-3       2-0   -1.612743304280462e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-3       2-1    1.178141724896658e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-3       2-2    4.297719780578690e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       2-3       2-3   -3.276060630732563e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-3       2-0   -1.948199092010796e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-3       2-1    2.249145748800627e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-3       2-2    9.384490959350733e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       2-3       2-3   -6.246086957622560e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-0       2-0    5.603508807100631e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-0       2-1    3.011990983999247e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-0       2-2    7.865604329382983e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-0       2-3    1.168101584347742e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-0       2-0   -4.813650081948034e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-0       2-1   -2.036304517101396e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-0       2-2   -6.820273066847479e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-0       2-3   -1.612743304277353e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-0       2-0    6.773888124661799e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-0       2-1    1.013387778081764e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-0       2-2    9.343503198756684e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-0       2-3   -1.911127836817447e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-0       2-0    5.832026133806252e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-0       2-1    3.369402767125522e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-0       2-2    8.263458600488613e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-0       2-3    1.954150468990917e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-1       2-0    3.011990983999370e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-1       2-1    1.667065778132682e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-1       2-2    4.422417099644329e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-1       2-3   -1.342531342686826e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-1       2-0   -2.036304517100828e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-1       2-1   -1.441006732780512e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-1       2-2   -3.441438905759958e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-1       2-3    1.178141724896658e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-1       2-0    1.013387778081796e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-1       2-1    1.992027385937283e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-1       2-2    1.474704180465490e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-1       2-3   -1.558268003342519e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-1       2-0    3.369402767126408e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-1       2-1    1.745906781681284e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-1       2-2    5.438852922658194e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-1       2-3   -1.427503855643367e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-2       2-0    7.865604329382985e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-2       2-1    4.422417099644103e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-2       2-2    2.737263975189050e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-2       2-3   -5.515976450460966e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-2       2-0   -6.820273066847479e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-2       2-1   -3.441438905759958e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-2       2-2   -2.386061760883756e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-2       2-3    4.297719780596454e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-2       2-0    9.343503198756682e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-2       2-1    1.474704180465465e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-2       2-2    3.218845468010853e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-2       2-3   -1.838243804837722e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-2       2-0    8.263458600488615e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-2       2-1    5.438852922644233e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-2       2-2    2.891016328506175e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-2       2-3   -6.748115547790645e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-3       2-0    1.168101584347759e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-3       2-1   -1.342531342686826e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-3       2-2   -5.515976450459634e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       2-3       2-3    3.728313871997137e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-3       2-0   -1.612743304278241e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-3       2-1    1.178141724896658e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-3       2-2    4.297719780596454e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       2-3       2-3   -3.276060630732563e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-3       2-0   -1.911127836816059e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-3       2-1   -1.558268003342519e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-3       2-2   -1.838243804837764e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       2-3       2-3    4.316313990912428e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-3       2-0    1.954150468993915e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-3       2-1   -1.427503855643367e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-3       2-2   -6.748115547790645e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       2-3       2-3    3.969485121903584e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-0       2-0    5.803069401295358e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-0       2-1    8.783173312644657e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-0       2-2    8.136045835280145e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-0       2-3    1.000086050382360e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-0       2-0   -9.386389032235364e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-0       2-1   -5.127844970083558e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-0       2-2   -1.317634215388623e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-0       2-3   -1.948199091997196e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-0       2-0    5.832026133806252e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-0       2-1    3.369402767125523e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-0       2-2    8.263458600488615e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-0       2-3    1.954150468983784e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-0       2-0    1.125655378841849e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-0       2-1    1.377880480183852e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-0       2-2    1.560964887384688e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-0       2-3    4.591055031437710e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-1       2-0    8.783173312644215e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-1       2-1    1.725072621922416e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-1       2-2    1.258453594401986e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-1       2-3   -1.386618534247727e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-1       2-0   -5.127844970083281e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-1       2-1   -2.792591155995217e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-1       2-2   -7.525499252661140e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-1       2-3    2.249145748800627e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-1       2-0    3.369402767126412e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-1       2-1    1.745906781681284e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-1       2-2    5.438852922658249e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-1       2-3   -1.427503855643367e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-1       2-0    1.377880480183849e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-1       2-1    3.321956681549609e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-1       2-2    2.007033777684852e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-1       2-3   -2.621999950836959e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-2       2-0    8.136045835280143e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-2       2-1    1.258453594401810e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-2       2-2    2.829543743995077e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-2       2-3   -1.547260227536812e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-2       2-0   -1.317634215388623e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-2       2-1   -7.525499252661029e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-2       2-2   -4.585562739975723e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-2       2-3    9.384490959358061e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-2       2-0    8.263458600488613e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-2       2-1    5.438852922658444e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-2       2-2    2.891016328506175e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-2       2-3   -6.748115547790867e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-2       2-0    1.560964887384688e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-2       2-1    2.007033777684791e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-2       2-2    5.394290296641206e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-2       2-3   -2.501691593554309e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-3       2-0    1.000086050382014e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-3       2-1   -1.386618534247728e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-3       2-2   -1.547260227536870e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       2-3       2-3    3.850284070723646e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-3       2-0   -1.948199092005522e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-3       2-1    2.249145748800627e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-3       2-2    9.384490959358727e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       2-3       2-3   -6.246086957622561e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-3       2-0    1.954150468993887e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-3       2-1   -1.427503855643367e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-3       2-2   -6.748115547791089e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       2-3       2-3    3.969485121903584e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-3       2-0    4.591055031524099e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-3       2-1   -2.621999950836958e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-3       2-2   -2.501691593554423e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       2-3       2-3    7.268588589507098e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-0       3-0    4.868192905220942e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-0       3-1   -4.978941720683451e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-0       3-2    6.786900972983623e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-0       3-3    9.090231107370371e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-0       3-0   -9.365712741203152e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-0       3-1   -1.224187107454375e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-0       3-2   -1.319230888499080e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-0       3-3   -5.024118596348992e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-0       3-0    1.704806367246882e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-0       3-1   -5.618522809324682e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-0       3-2    2.391797155811619e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-0       3-3    1.823456440428963e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-0       3-0    1.770608736627365e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-0       3-1   -1.775649667517792e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-0       3-2    2.478890866452099e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-0       3-3    2.429558215550536e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-1       3-0   -4.978941720684014e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-1       3-1    1.444686338075696e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-1       3-2   -7.506729201028899e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-1       3-3   -1.150521248265149e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-1       3-0   -1.224187107454375e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-1       3-1   -2.798483263165502e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-1       3-2   -2.968168669821192e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-1       3-3    2.266151169866135e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-1       3-0   -5.618522809322488e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-1       3-1    5.080483535515440e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-1       3-2   -7.248843993927223e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-1       3-3   -4.087806245312244e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-1       3-0   -1.775649667518070e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-1       3-1    5.269197572866622e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-1       3-2   -2.528297801328605e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-1       3-3   -4.224319830024649e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-2       3-0    6.786900972983625e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-2       3-1   -7.506729201029664e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-2       3-2    2.363302657904221e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-2       3-3    9.678660751827259e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-2       3-0   -1.319230888499080e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-2       3-1   -2.968168669823412e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-2       3-2   -4.620440742267817e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-2       3-3    4.280831844716865e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-2       3-0    2.391797155811619e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-2       3-1   -7.248843993930554e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-2       3-2    8.358364188694710e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-2       3-3    8.858804397942482e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-2       3-0    2.478890866452099e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-2       3-1   -2.528297801328154e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-2       3-2    8.651367387525942e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-2       3-3    3.174495389033239e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-3       3-0    9.090231107375911e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-3       3-1   -1.150521248265149e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-3       3-2    9.678660751827536e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       3-3       3-3    3.203183089524514e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-3       3-0   -5.024118596340110e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-3       3-1    2.266151169866135e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-3       3-2    4.280831844715977e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       3-3       3-3   -6.317271249999760e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-3       3-0    1.823456440431183e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-3       3-1   -4.087806245312244e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-3       3-2    8.858804397936376e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       3-3       3-3    1.138996504616417e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-3       3-0    2.429558215549148e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-3       3-1   -4.224319830024649e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-3       3-2    3.174495389033230e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       3-3       3-3    1.176485284577908e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-0       3-0   -9.365712741203152e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-0       3-1   -1.224187107454542e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-0       3-2   -1.319230888499080e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-0       3-3   -5.024118596344551e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-0       3-0    1.440631927265614e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-0       3-1   -1.484671096951944e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-0       3-2    2.008677179921740e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-0       3-3    2.708366727979580e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-0       3-0   -1.307771396212538e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-0       3-1   -1.594355395313801e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-0       3-2   -1.842113164955113e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-0       3-3   -7.003157820455880e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-0       3-0   -2.843362947214892e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-0       3-1    9.938705110478789e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-0       3-2   -3.989390222975097e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-0       3-3   -3.088119259200368e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-1       3-0   -1.224187107454486e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-1       3-1   -2.798483263165502e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-1       3-2   -2.968168669822080e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-1       3-3    2.266151169866135e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-1       3-0   -1.484671096952208e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-1       3-1    4.275572290648780e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-1       3-2   -2.240451139876836e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-1       3-3   -3.405711808682705e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-1       3-0   -1.594355395313884e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-1       3-1   -3.907659430207085e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-1       3-2   -3.984348402467155e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-1       3-3    3.164374490737676e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-1       3-0    9.938705110478861e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-1       3-1   -8.473812976653097e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-1       3-2    1.297154088092951e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-1       3-3    6.818783893709981e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-2       3-0   -1.319230888499080e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-2       3-1   -2.968168669824300e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-2       3-2   -4.620440742267817e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-2       3-3    4.280831844719529e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-2       3-0    2.008677179921741e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-2       3-1   -2.240451139876829e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-2       3-2    6.995058956790282e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-2       3-3    2.890478137018426e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-2       3-0   -1.842113164955113e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-2       3-1   -3.984348402467162e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-2       3-2   -6.451787876375055e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-2       3-3    5.780791167357080e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-2       3-0   -3.989390222975097e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-2       3-1    1.297154088092466e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-2       3-2   -1.394180390478837e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-2       3-3   -1.592907766039570e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-3       3-0   -5.024118596335669e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-3       3-1    2.266151169866135e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-3       3-2    4.280831844718641e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       3-3       3-3   -6.317271249999762e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-3       3-0    2.708366727981943e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-3       3-1   -3.405711808682705e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-3       3-2    2.890478137017427e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       3-3       3-3    9.482150203619428e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-3       3-0   -7.003157820471423e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-3       3-1    3.164374490737676e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-3       3-2    5.780791167356414e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       3-3       3-3   -8.821185713791088e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-3       3-0   -3.088119259201964e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-3       3-1    6.818783893709981e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-3       3-2   -1.592907766039126e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       3-3       3-3   -1.899960295526905e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-0       3-0    1.704806367246882e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-0       3-1   -5.618522809324127e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-0       3-2    2.391797155811619e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-0       3-3    1.823456440431183e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-0       3-0   -1.307771396212538e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-0       3-1   -1.594355395314689e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-0       3-2   -1.842113164955113e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-0       3-3   -7.003157820446998e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-0       3-0    2.458529777497663e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-0       3-1   -2.444421453804678e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-0       3-2    3.430117604968551e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-0       3-3    4.512576863992561e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-0       3-0    1.583476247310764e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-0       3-1    1.761193592760784e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-0       3-2    2.230426808342835e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-0       3-3    8.516818356716326e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-1       3-0   -5.618522809321377e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-1       3-1    5.080483535515441e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-1       3-2   -7.248843993932774e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-1       3-3   -4.087806245312245e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-1       3-0   -1.594355395312108e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-1       3-1   -3.907659430207085e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-1       3-2   -3.984348402467155e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-1       3-3    3.164374490737676e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-1       3-0   -2.444421453804789e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-1       3-1    7.299618975045829e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-1       3-2   -3.678331749636781e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-1       3-3   -5.820642824113400e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-1       3-0    1.761193592759313e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-1       3-1    4.731415219420423e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-1       3-2    4.578930820003289e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-1       3-3   -3.831309181077405e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-2       3-0    2.391797155811619e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-2       3-1   -7.248843993932774e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-2       3-2    8.358364188694708e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-2       3-3    8.858804397942482e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-2       3-0   -1.842113164955113e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-2       3-1   -3.984348402467162e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-2       3-2   -6.451787876375055e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-2       3-3    5.780791167357080e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-2       3-0    3.430117604968551e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-2       3-1   -3.678331749636101e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-2       3-2    1.194950806674837e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-2       3-3    4.742933904569140e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-2       3-0    2.230426808342835e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-2       3-1    4.578930820006787e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-2       3-2    7.811707587205272e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-2       3-3   -6.692123590079924e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-3       3-0    1.823456440432293e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-3       3-1   -4.087806245312245e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-3       3-2    8.858804397945258e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       3-3       3-3    1.138996504616417e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-3       3-0   -7.003157820453659e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-3       3-1    3.164374490737676e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-3       3-2    5.780791167356414e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       3-3       3-3   -8.821185713791088e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-3       3-0    4.512576863993046e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-3       3-1   -5.820642824113400e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-3       3-2    4.742933904569806e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       3-3       3-3    1.620727607236622e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-3       3-0    8.516818356711885e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-3       3-1   -3.831309181077405e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-3       3-2   -6.692123590081256e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       3-3       3-3    1.068029552075303e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-0       3-0    1.770608736627365e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-0       3-1   -1.775649667517128e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-0       3-2    2.478890866452101e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-0       3-3    2.429558215571271e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-0       3-0   -2.843362947214892e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-0       3-1    9.938705110477614e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-0       3-2   -3.989390222975098e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-0       3-3   -3.088119259201782e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-0       3-0    1.583476247310764e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-0       3-1    1.761193592760785e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-0       3-2    2.230426808342835e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-0       3-3    8.516818356716291e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-0       3-0    3.864160520732090e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-0       3-1   -3.293648233935837e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-0       3-2    5.399440002614476e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-0       3-3    6.422385887872679e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-1       3-0   -1.775649667517232e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-1       3-1    5.269197572866624e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-1       3-2   -2.528297801330567e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-1       3-3   -4.224319830024655e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-1       3-0    9.938705110475394e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-1       3-1   -8.473812976653099e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-1       3-2    1.297154088092067e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-1       3-3    6.818783893709982e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-1       3-0    1.761193592759314e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-1       3-1    4.731415219420423e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-1       3-2    4.578930820003303e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-1       3-3   -3.831309181077405e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-1       3-0   -3.293648233935774e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-1       3-1    1.148466072039579e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-1       3-2   -4.890999841216086e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-1       3-3   -9.180506854963873e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-2       3-0    2.478890866452101e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-2       3-1   -2.528297801331564e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-2       3-2    8.651367387525950e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-2       3-3    3.174495389037437e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-2       3-0   -3.989390222975098e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-2       3-1    1.297154088092247e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-2       3-2   -1.394180390478837e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-2       3-3   -1.592907766039153e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-2       3-0    2.230426808342835e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-2       3-1    4.578930820006794e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-2       3-2    7.811707587205272e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-2       3-3   -6.692123590072874e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-2       3-0    5.399440002614476e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-2       3-1   -4.890999841216100e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-2       3-2    1.882624769702501e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-2       3-3    6.281346739047594e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-3       3-0    2.429558215571098e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-3       3-1   -4.224319830024655e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-3       3-2    3.174495389036743e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       3-3       3-3    1.176485284577909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-3       3-0   -3.088119259200880e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-3       3-1    6.818783893709981e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-3       3-2   -1.592907766038279e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       3-3       3-3   -1.899960295526905e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-3       3-0    8.516818356711850e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-3       3-1   -3.831309181077405e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-3       3-2   -6.692123590081311e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       3-3       3-3    1.068029552075303e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-3       3-0    6.422385887873226e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-3       3-1   -9.180506854963873e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-3       3-2    6.281346739051036e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       3-3       3-3    2.556761666809046e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-0       4-0   -3.465293304441793e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-0       4-1   -4.765379121600360e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-0       4-2   -3.529021706453871e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-0       4-3    1.503826760371886e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-0       4-0   -1.098527899698855e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-0       4-1   -7.309281477041767e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-0       4-2   -1.527006930018747e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-0       4-3   -4.800510140058125e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-0       4-0    7.341476786867059e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-0       4-1   -6.689517381855909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-0       4-2    1.202650052050388e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-0       4-3    2.688733192835091e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-0       4-0    1.702575395138278e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-0       4-1   -7.538069686394482e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-0       4-2    4.737776666463243e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-0       4-3    1.144973230225925e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-1       4-0   -4.765379121600360e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-1       4-1   -1.020638378788489e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-1       4-2   -7.003383936804995e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-1       4-3    5.835790438901960e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-1       4-0   -7.309281477041762e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-1       4-1   -3.254905437234671e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-1       4-2   -1.250697603295910e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-1       4-3    2.581653507110617e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-1       4-0   -6.689517381855909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-1       4-1    2.182750253605945e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-1       4-2   -9.814266968538595e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-1       4-3   -2.050091228720545e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-1       4-0   -7.538069686394482e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-1       4-1    5.596521957627279e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-1       4-2   -1.079583482943085e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-1       4-3   -8.781300996130030e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-2       4-0   -3.529021706453871e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-2       4-1   -7.003383936804995e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-2       4-2   -1.630367089577443e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-2       4-3    8.822740556184289e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-2       4-0   -1.527006930018747e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-2       4-1   -1.250697603295905e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-2       4-2   -5.313175139925755e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-2       4-3    1.472641299155634e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-2       4-0    1.202650052050389e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-2       4-1   -9.814266968538595e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-2       4-2    3.611156291969887e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-2       4-3    1.237874671141335e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-2       4-0    4.737776666463243e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-2       4-1   -1.079583482943084e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-2       4-2    1.060116058938348e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-2       4-3    1.331881796054531e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-3       4-0    1.503826760371887e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-3       4-1    5.835790438901964e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-3       4-2    8.822740556184289e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       4-3       4-3   -2.101186224217623e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-3       4-0   -4.800510140058116e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-3       4-1    2.581653507110617e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-3       4-2    1.472641299155632e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       4-3       4-3   -7.188057725291665e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-3       4-0    2.688733192835088e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-3       4-1   -2.050091228720545e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-3       4-2    1.237874671141335e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       4-3       4-3    5.030772821114048e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-3       4-0    1.144973230225932e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-3       4-1   -8.781300996130031e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-3       4-2    1.331881796054531e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       4-3       4-3    1.698737086343529e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-0       4-0   -1.098527899698855e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-0       4-1   -7.309281477041776e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-0       4-2   -1.527006930018747e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-0       4-3   -4.800510140058125e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-0       4-0   -1.020545621396667e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-0       4-1   -1.435757900674959e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-0       4-2   -1.031747846686775e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-0       4-3    4.564791190291607e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-0       4-0   -1.538207264802467e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-0       4-1   -1.704320154237142e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-0       4-2   -2.136912507867336e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-0       4-3   -6.209952722984738e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-0       4-0   -1.215041892865140e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-0       4-1    1.158705074860311e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-0       4-2   -2.003442887501267e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-0       4-3   -4.621414183516382e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-1       4-0   -7.309281477041771e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-1       4-1   -3.254905437234671e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-1       4-2   -1.250697603295912e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-1       4-3    2.581653507110618e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-1       4-0   -1.435757900674959e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-1       4-1   -3.007128623086227e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-1       4-2   -2.110412443370667e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-1       4-3    1.707472803201358e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-1       4-0   -1.704320154237141e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-1       4-1   -4.558394421831071e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-1       4-2   -2.798566044656384e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-1       4-3    3.613635096754052e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-1       4-0    1.158705074860311e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-1       4-1   -3.611799055348516e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-1       4-2    1.699880529137676e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-1       4-3    3.414092875086007e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-2       4-0   -1.527006930018747e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-2       4-1   -1.250697603295909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-2       4-2   -5.313175139925755e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-2       4-3    1.472641299155637e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-2       4-0   -1.031747846686776e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-2       4-1   -2.110412443370667e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-2       4-2   -4.804880305456107e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-2       4-3    2.659174931460785e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-2       4-0   -2.136912507867336e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-2       4-1   -2.798566044656384e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-2       4-2   -7.442051247128744e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-2       4-3    3.432577201636957e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-2       4-0   -2.003442887501267e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-2       4-1    1.699880529137676e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-2       4-2   -5.976050643525406e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-2       4-3   -2.143930611130103e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-3       4-0   -4.800510140058125e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-3       4-1    2.581653507110618e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-3       4-2    1.472641299155632e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       4-3       4-3   -7.188057725291665e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-3       4-0    4.564791190291609e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-3       4-1    1.707472803201359e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-3       4-2    2.659174931460785e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       4-3       4-3   -6.190022673150614e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-3       4-0   -6.209952722984756e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-3       4-1    3.613635096754052e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-3       4-2    3.432577201636955e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       4-3       4-3   -1.006844808910578e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-3       4-0   -4.621414183516380e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-3       4-1    3.414092875086007e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-3       4-2   -2.143930611130103e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       4-3       4-3   -8.332465158669056e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-0       4-0    7.341476786867058e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-0       4-1   -6.689517381855911e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-0       4-2    1.202650052050389e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-0       4-3    2.688733192835089e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-0       4-0   -1.538207264802467e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-0       4-1   -1.704320154237142e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-0       4-2   -2.136912507867336e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-0       4-3   -6.209952722984738e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-0       4-0   -1.512771642312835e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-0       4-1   -2.388041710644127e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-0       4-2   -1.461174223657409e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-0       4-3    7.701854407423677e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-0       4-0    1.858851468877211e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-0       4-1    2.431590476443080e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-0       4-2    2.582024897568699e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-0       4-3    7.048831308827467e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-1       4-0   -6.689517381855909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-1       4-1    2.182750253605945e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-1       4-2   -9.814266968538593e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-1       4-3   -2.050091228720545e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-1       4-0   -1.704320154237142e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-1       4-1   -4.558394421831071e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-1       4-2   -2.798566044656387e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-1       4-3    3.613635096754052e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-1       4-0   -2.388041710644127e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-1       4-1   -4.457487436965645e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-1       4-2   -3.510306119103466e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-1       4-3    2.413531622376867e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-1       4-0    2.431590476443080e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-1       4-1    5.509516266105675e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-1       4-2    3.984598866395940e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-1       4-3   -4.367436746940858e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-2       4-0    1.202650052050390e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-2       4-1   -9.814266968538593e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-2       4-2    3.611156291969888e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-2       4-3    1.237874671141335e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-2       4-0   -2.136912507867336e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-2       4-1   -2.798566044656384e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-2       4-2   -7.442051247128744e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-2       4-3    3.432577201636957e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-2       4-0   -1.461174223657409e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-2       4-1   -3.510306119103466e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-2       4-2   -7.110023587662028e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-2       4-3    4.423849446306932e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-2       4-0    2.582024897568699e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-2       4-1    3.984598866395942e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-2       4-2    8.996475260148989e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-2       4-3   -4.967828328839572e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-3       4-0    2.688733192835089e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-3       4-1   -2.050091228720545e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-3       4-2    1.237874671141335e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       4-3       4-3    5.030772821114048e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-3       4-0   -6.209952722984738e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-3       4-1    3.613635096754052e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-3       4-2    3.432577201636955e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       4-3       4-3   -1.006844808910578e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-3       4-0    7.701854407423671e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-3       4-1    2.413531622376866e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-3       4-2    4.423849446306932e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       4-3       4-3   -9.110711487892390e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-3       4-0    7.048831308827472e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-3       4-1   -4.367436746940858e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-3       4-2   -4.967828328839571e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       4-3       4-3    1.217271776215927e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-0       4-0    1.702575395138262e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-0       4-1   -7.538069686394484e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-0       4-2    4.737776666463232e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-0       4-3    1.144973230225942e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-0       4-0   -1.215041892865139e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-0       4-1    1.158705074860310e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-0       4-2   -2.003442887501266e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-0       4-3   -4.621414183516380e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-0       4-0    1.858851468877211e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-0       4-1    2.431590476443077e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-0       4-2    2.582024897568699e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-0       4-3    7.048831308827471e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-0       4-0   -1.354121964237308e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-0       4-1   -3.263598984456568e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-0       4-2   -1.005835211001091e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-0       4-3    1.083704108835118e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-1       4-0   -7.538069686394485e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-1       4-1    5.596521957627253e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-1       4-2   -1.079583482943085e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-1       4-3   -8.781300996130031e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-1       4-0    1.158705074860310e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-1       4-1   -3.611799055348513e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-1       4-2    1.699880529137676e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-1       4-3    3.414092875086007e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-1       4-0    2.431590476443077e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-1       4-1    5.509516266105673e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-1       4-2    3.984598866395940e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-1       4-3   -4.367436746940858e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-1       4-0   -3.263598984456568e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-1       4-1   -3.980706906346838e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-1       4-2   -4.795904614867999e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-1       4-3    1.626796900492233e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-2       4-0    4.737776666463230e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-2       4-1   -1.079583482943085e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-2       4-2    1.060116058938346e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-2       4-3    1.331881796054532e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-2       4-0   -2.003442887501266e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-2       4-1    1.699880529137676e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-2       4-2   -5.976050643525406e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-2       4-3   -2.143930611130103e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-2       4-0    2.582024897568699e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-2       4-1    3.984598866395939e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-2       4-2    8.996475260148989e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-2       4-3   -4.967828328839571e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-2       4-0   -1.005835211001092e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-2       4-1   -4.795904614867999e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-2       4-2   -6.276124942444375e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-2       4-3    6.044568179828537e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-3       4-0    1.144973230225940e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-3       4-1   -8.781300996130041e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-3       4-2    1.331881796054532e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       4-3       4-3    1.698737086343531e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-3       4-0   -4.621414183516381e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-3       4-1    3.414092875086005e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-3       4-2   -2.143930611130103e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       4-3       4-3   -8.332465158669056e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-3       4-0    7.048831308827472e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-3       4-1   -4.367436746940858e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-3       4-2   -4.967828328839577e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       4-3       4-3    1.217271776215927e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-3       4-0    1.083704108835119e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-3       4-1    1.626796900492230e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-3       4-2    6.044568179828536e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       4-3       4-3   -7.802680766640719e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-0       5-0    1.507668687586465e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-0       5-1   -1.710530491573905e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-0       5-2    2.001414518666849e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-0       5-3   -1.747146290841908e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-0       5-0   -1.859042103140767e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-0       5-1   -8.147256170494338e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-0       5-2   -3.047666147176431e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-0       5-3    2.382666866030106e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-0       5-0    2.280934998243985e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-0       5-1   -2.458547511818383e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-0       5-2    3.092107961539629e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-0       5-3   -5.589892621061035e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-0       5-0    2.265118141249791e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-0       5-1   -2.503096611021207e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-0       5-2    3.023104191460079e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-0       5-3    2.368125507219194e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-1       5-0   -1.710530491574014e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-1       5-1    4.333379554149221e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-1       5-2   -2.379320423558250e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-1       5-3   -3.171381542716606e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-1       5-0   -8.147256173061728e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-1       5-1   -6.162826339642987e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-1       5-2   -4.199212446306822e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-1       5-3    6.172501648478967e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-1       5-0   -2.458547511818439e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-1       5-1    6.646650115613724e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-1       5-2   -3.366211207077406e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-1       5-3   -5.048191559118909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-1       5-0   -2.503096611021212e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-1       5-1    6.533337970260241e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-1       5-2   -7.518332388126278e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-1       5-3   -4.831585065598937e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-2       5-0    2.001414518666849e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-2       5-1   -2.379320423558218e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-2       5-2    6.773527458671185e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-2       5-3    2.862765998960872e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-2       5-0   -3.047666147176431e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-2       5-1   -4.199212446318080e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-2       5-2   -1.151922961730538e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-2       5-3    9.596973145066384e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-2       5-0    3.092107961539630e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-2       5-1   -3.366211207077851e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-2       5-2    1.059799670419278e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-2       5-3    3.991479357456300e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-2       5-0    3.023104191460080e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-2       5-1   -7.518332388133217e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-2       5-2    1.026956618491721e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-2       5-3    1.355216271128040e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-3       5-0   -1.747146290862747e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-3       5-1   -3.171381542716606e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-3       5-2    2.862765998960872e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-0       5-3       5-3    8.770110690627037e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-3       5-0    2.382666866029906e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-3       5-1    6.172501648478966e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-3       5-2    9.596973145048482e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-1       5-3       5-3   -1.746327086687586e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-3       5-0   -5.589892621060688e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-3       5-1   -5.048191559118909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-3       5-2    3.991479357454406e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-2       5-3       5-3    1.400334480378907e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-3       5-0    2.368125507219709e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-3       5-1   -4.831585065598937e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-3       5-2    1.355216271128868e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-0       1-3       5-3       5-3    1.338228166136394e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-0       5-0   -1.859042103140767e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-0       5-1   -8.147256169473019e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-0       5-2   -3.047666147176431e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-0       5-3    2.382666866030470e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-0       5-0    4.449255298082021e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-0       5-1   -5.194854944889376e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-0       5-2    5.906817744616427e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-0       5-3   -5.958243606481433e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-0       5-0   -2.355782253660207e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-0       5-1   -6.809601998762390e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-0       5-2   -3.932926138304508e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-0       5-3    2.804161167427686e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-0       5-0   -3.781471303148679e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-0       5-1    4.241112011927422e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-0       5-2   -5.126075709012605e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-0       5-3    9.666398188229106e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-1       5-0   -8.147256172040410e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-1       5-1   -6.162826339642987e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-1       5-2   -4.199212446307533e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-1       5-3    6.172501648478967e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-1       5-0   -5.194854944889229e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-1       5-1    1.278883708728538e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-1       5-2   -7.215242889997486e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-1       5-3   -9.360612725811333e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-1       5-0   -6.809601998761835e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-1       5-1   -7.909784333882194e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-1       5-2   -1.455382595394777e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-1       5-3    8.099182913720057e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-1       5-0    4.241112011927207e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-1       5-1   -1.101892374439765e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-1       5-2    5.807075698186719e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-1       5-3    8.368095909637063e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-2       5-0   -3.047666147176431e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-2       5-1   -4.199212446318636e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-2       5-2   -1.151922961730538e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-2       5-3    9.596973145069160e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-2       5-0    5.906817744616427e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-2       5-1   -7.215242889997539e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-2       5-2    1.999152991872992e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-2       5-3    8.668872220303381e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-2       5-0   -3.932926138304508e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-2       5-1   -1.455382595395110e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-2       5-2   -1.498543080998993e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-2       5-3    2.300216498124424e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-2       5-0   -5.126075709012605e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-2       5-1    5.807075698187440e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-2       5-2   -1.756850691569377e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-2       5-3   -6.885350649264266e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-3       5-0    2.382666866030192e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-3       5-1    6.172501648478966e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-3       5-2    9.596973145051535e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-0       5-3       5-3   -1.746327086687586e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-3       5-0   -5.958243606482127e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-3       5-1   -9.360612725811333e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-3       5-2    8.668872220303381e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-1       5-3       5-3    2.588538626551485e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-3       5-0    2.804161167427755e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-3       5-1    8.099182913720059e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-3       5-2    2.300216498124424e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-2       5-3       5-3   -2.294660883494409e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-3       5-0    9.666398188248534e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-3       5-1    8.368095909637063e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-3       5-2   -6.885350649264890e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-1       1-3       5-3       5-3   -2.321162618495589e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-0       5-0    2.280934998243985e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-0       5-1   -2.458547511818438e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-0       5-2    3.092107961539630e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-0       5-3   -5.589892621055484e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-0       5-0   -2.355782253660207e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-0       5-1   -6.809601998761818e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-0       5-2   -3.932926138304508e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-0       5-3    2.804161167429837e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-0       5-0    7.212719601569674e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-0       5-1   -8.683774054916934e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-0       5-2    9.583114801324422e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-0       5-3   -1.117621167649255e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-0       5-0    2.709994426052127e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-0       5-1    1.436438987852199e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-0       5-2    4.566042485010723e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-0       5-3   -2.822357573080483e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-1       5-0   -2.458547511818495e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-1       5-1    6.646650115613725e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-1       5-2   -3.366211207078295e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-1       5-3   -5.048191559118909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-1       5-0   -6.809601998759597e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-1       5-1   -7.909784333882195e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-1       5-2   -1.455382595393875e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-1       5-3    8.099182913720055e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-1       5-0   -8.683774054916753e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-1       5-1    2.074271838525672e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-1       5-2   -1.204100221009869e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-1       5-3   -1.520342452181451e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-1       5-0    1.436438987852220e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-1       5-1    9.158124462552678e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-1       5-2    2.527106980328493e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-1       5-3   -9.480774057687130e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-2       5-0    3.092107961539630e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-2       5-1   -3.366211207078295e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-2       5-2    1.059799670419278e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-2       5-3    3.991479357456300e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-2       5-0   -3.932926138304508e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-2       5-1   -1.455382595393972e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-2       5-2   -1.498543080998993e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-2       5-3    2.300216498124202e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-2       5-0    9.583114801324422e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-2       5-1   -1.204100221009872e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-2       5-2    3.244891859863088e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-2       5-3    1.444382293238956e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-2       5-0    4.566042485010723e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-2       5-1    2.527106980328216e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-2       5-2    1.746782417121141e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-2       5-3   -3.609105900829218e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-3       5-0   -5.589892621032932e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-3       5-1   -5.048191559118909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-3       5-2    3.991479357456626e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-0       5-3       5-3    1.400334480378907e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-3       5-0    2.804161167428727e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-3       5-1    8.099182913720057e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-3       5-2    2.300216498124313e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-1       5-3       5-3   -2.294660883494409e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-3       5-0   -1.117621167653141e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-3       5-1   -1.520342452181451e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-3       5-2    1.444382293238958e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-2       5-3       5-3    4.204657368848392e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-3       5-0   -2.822357573078124e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-3       5-1   -9.480774057687130e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-3       5-2   -3.609105900831105e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-2       1-3       5-3       5-3    2.688050490233316e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-0       5-0    2.265118141249790e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-0       5-1   -2.503096611027327e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-0       5-2    3.023104191460078e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-0       5-3    2.368125507221536e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-0       5-0   -3.781471303148679e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-0       5-1    4.241112011927193e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-0       5-2   -5.126075709012605e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-0       5-3    9.666398188276290e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-0       5-0    2.709994426052127e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-0       5-1    1.436438987852254e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-0       5-2    4.566042485010724e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-0       5-3   -2.822357573082687e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-0       5-0    9.704767424135522e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-0       5-1   -1.185512729852781e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-0       5-2    1.292845924814172e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-0       5-3   -1.706330798085052e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-1       5-0   -2.503096611027333e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-1       5-1    6.533337970260237e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-1       5-2   -7.518332388141815e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-1       5-3   -4.831585065598937e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-1       5-0    4.241112011926978e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-1       5-1   -1.101892374439765e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-1       5-2    5.807075698186663e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-1       5-3    8.368095909637065e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-1       5-0    1.436438987852275e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-1       5-1    9.158124462552678e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-1       5-2    2.527106980328715e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-1       5-3   -9.480774057687130e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-1       5-0   -1.185512729852735e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-1       5-1    2.795797167323528e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-1       5-2   -1.640782063090800e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-1       5-3   -2.058964955713594e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-2       5-0    3.023104191460078e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-2       5-1   -7.518332388151442e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-2       5-2    1.026956618491721e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-2       5-3    1.355216271128772e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-2       5-0   -5.126075709012605e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-2       5-1    5.807075698184665e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-2       5-2   -1.756850691569377e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-2       5-3   -6.885350649263822e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-2       5-0    4.566042485010724e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-2       5-1    2.527106980328882e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-2       5-2    1.746782417121141e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-2       5-3   -3.609105900829218e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-2       5-0    1.292845924814172e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-2       5-1   -1.640782063090803e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-2       5-2    4.384702523137884e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-2       5-3    1.964778260872771e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-3       5-0    2.368125507220712e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-3       5-1   -4.831585065598936e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-3       5-2    1.355216271128491e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-0       5-3       5-3    1.338228166136393e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-3       5-0    9.666398188251310e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-3       5-1    8.368095909637063e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-3       5-2   -6.885350649264668e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-1       5-3       5-3   -2.321162618495589e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-3       5-0   -2.822357573080327e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-3       5-1   -9.480774057687130e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-3       5-2   -3.609105900831078e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-2       5-3       5-3    2.688050490233316e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-3       5-0   -1.706330798085052e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-3       5-1   -2.058964955713594e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-3       5-2    1.964778260872749e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "1-3       1-3       5-3       5-3    5.696461764951018e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-0       3-0    8.593696846210803e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-0       3-1   -3.712726902540870e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-0       3-2    1.190043708278496e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-0       3-3    1.908391885239490e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-0       3-0   -7.631281043109200e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-0       3-1   -7.068445733541373e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-0       3-2   -1.088457603893902e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-0       3-3    2.774554798583352e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-0       3-0    1.191120674271124e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-0       3-1   -5.307289985962373e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-0       3-2    1.649362138566599e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-0       3-3    2.602905336256572e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-0       3-0   -7.404741688721612e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-0       3-1    1.621239212646677e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-0       3-2   -1.170105244613001e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-0       3-3    3.962265087723473e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-1       3-0   -3.712726902540869e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-1       3-1    2.538921679842137e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-1       3-2   -5.378410847611535e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-1       3-3   -1.999330371351610e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-1       3-0   -7.068445733541373e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-1       3-1   -2.305496487916245e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-1       3-2   -1.474140924754372e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-1       3-3    1.902792930897357e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-1       3-0   -5.307289985962373e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-1       3-1    3.518925345829124e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-1       3-2   -7.675116399889866e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-1       3-3   -2.770816186351187e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-1       3-0    1.621239212646677e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-1       3-1   -2.253383717000880e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-1       3-2    1.607220219523157e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-1       3-3    1.995546275556981e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-2       3-0    1.190043708278496e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-2       3-1   -5.378410847611538e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-2       3-2    4.127630842205885e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-2       3-3    6.761312877254805e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-2       3-0   -1.088457603893903e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-2       3-1   -1.474140924754372e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-2       3-2   -3.856762195665804e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-2       3-3    2.235085579994988e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-2       3-0    1.649362138566599e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-2       3-1   -7.675116399889964e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-2       3-2    5.720591551481931e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-2       3-3    9.633770415046608e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-2       3-0   -1.170105244613059e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-2       3-1    1.607220219523157e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-2       3-2   -3.776084548214787e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-2       3-3   -1.181958275632973e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-3       3-0    1.908391885239473e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-3       3-1   -1.999330371351610e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-3       3-2    6.761312877254847e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       3-3       3-3    5.560614523863995e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-3       3-0    2.774554798583351e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-3       3-1    1.902792930897351e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-3       3-2    2.235085579994988e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       3-3       3-3   -5.329251342734319e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-3       3-0    2.602905336257107e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-3       3-1   -2.770816186351187e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-3       3-2    9.633770415047049e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       3-3       3-3    7.706246876744716e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-3       3-0    3.962265087723473e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-3       3-1    1.995546275556982e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-3       3-2   -1.181958275632973e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       3-3       3-3   -5.204268186958902e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-0       3-0   -7.631281043109200e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-0       3-1   -7.068445733541371e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-0       3-2   -1.088457603893908e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-0       3-3    2.774554798583351e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-0       3-0    2.535293413426937e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-0       3-1   -1.119880425823585e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-0       3-2    3.510719573745098e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-0       3-3    5.642555573757433e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-0       3-0   -1.055278887179905e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-0       3-1   -1.271396645361700e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-0       3-2   -1.503450971924459e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-0       3-3    3.226726226843021e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-0       3-0   -1.999272217942325e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-0       3-1    9.292339866626427e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-0       3-2   -2.768228497911302e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-0       3-3   -4.409689592288538e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-1       3-0   -7.068445733541373e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-1       3-1   -2.305496487916267e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-1       3-2   -1.474140924754372e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-1       3-3    1.902792930897369e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-1       3-0   -1.119880425823596e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-1       3-1    7.490099515386965e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-1       3-2   -1.621606165251367e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-1       3-3   -5.897907679743459e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-1       3-0   -1.271396645361700e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-1       3-1   -3.188226517015873e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-1       3-2   -2.348516580790067e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-1       3-3    2.629669015511564e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-1       3-0    9.292339866626080e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-1       3-1   -5.906173289163096e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-1       3-2    1.343175849379120e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-1       3-3    4.650006550235300e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-2       3-0   -1.088457603893909e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-2       3-1   -1.474140924754372e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-2       3-2   -3.856762195665850e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-2       3-3    2.235085579994988e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-2       3-0    3.510719573745098e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-2       3-1   -1.621606165251356e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-2       3-2    1.217658384893180e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-2       3-3    2.037461977287205e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-2       3-0   -1.503450971924451e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-2       3-1   -2.348516580790067e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-2       3-2   -5.334141491260880e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-2       3-3    3.346440636973141e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-2       3-0   -2.768228497911302e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-2       3-1    1.343175849379110e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-2       3-2   -9.600844695299751e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-2       3-3   -1.684755883826707e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-3       3-0    2.774554798583351e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-3       3-1    1.902792930897363e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-3       3-2    2.235085579994988e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       3-3       3-3   -5.329251342734372e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-3       3-0    5.642555573755767e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-3       3-1   -5.897907679743457e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-3       3-2    2.037461977287310e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       3-3       3-3    1.640341208373102e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-3       3-0    3.226726226843021e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-3       3-1    2.629669015511541e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-3       3-2    3.346440636973141e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       3-3       3-3   -7.371885461052714e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-3       3-0   -4.409689592288538e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-3       3-1    4.650006550235300e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-3       3-2   -1.684755883826895e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       3-3       3-3   -1.293256617135361e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-0       3-0    1.191120674271124e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-0       3-1   -5.307289985962318e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-0       3-2    1.649362138566599e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-0       3-3    2.602905336256016e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-0       3-0   -1.055278887179905e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-0       3-1   -1.271396645361700e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-0       3-2   -1.503450971924471e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-0       3-3    3.226726226843021e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-0       3-0    4.114808070987604e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-0       3-1   -1.872493863169646e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-0       3-2    5.697661690107226e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-0       3-3    9.193625318139539e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-0       3-0    1.261566154344850e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-0       3-1    1.771078435193314e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-0       3-2    1.795899585901107e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-0       3-3   -3.144563100581599e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-1       3-0   -5.307289985962318e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-1       3-1    3.518925345829124e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-1       3-2   -7.675116399889645e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-1       3-3   -2.770816186351188e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-1       3-0   -1.271396645361700e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-1       3-1   -3.188226517015896e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-1       3-2   -2.348516580790067e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-1       3-3    2.629669015511570e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-1       3-0   -1.872493863169635e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-1       3-1    1.215612459756988e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-1       3-2   -2.710008937796015e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-1       3-3   -9.571299271674082e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-1       3-0    1.771078435193314e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-1       3-1    3.811566183504587e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-1       3-2    3.039539191280858e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-1       3-3   -3.142728259861854e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-2       3-0    1.649362138566599e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-2       3-1   -7.675116399889742e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-2       3-2    5.720591551481931e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-2       3-3    9.633770415046608e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-2       3-0   -1.503450971924462e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-2       3-1   -2.348516580790067e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-2       3-2   -5.334141491260886e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-2       3-3    3.346440636973141e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-2       3-0    5.697661690107226e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-2       3-1   -2.710008937796103e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-2       3-2    1.976123084210368e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-2       3-3    3.402773186359155e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-2       3-0    1.795899585901104e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-2       3-1    3.039539191280858e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-2       3-2    6.377867042855767e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-2       3-3   -4.140640161473054e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-3       3-0    2.602905336256279e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-3       3-1   -2.770816186351188e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-3       3-2    9.633770415046834e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       3-3       3-3    7.706246876744716e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-3       3-0    3.226726226843021e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-3       3-1    2.629669015511569e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-3       3-2    3.346440636973141e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       3-3       3-3   -7.371885461052724e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-3       3-0    9.193625318138421e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-3       3-1   -9.571299271674082e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-3       3-2    3.402773186359098e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       3-3       3-3    2.661977178174370e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-3       3-0   -3.144563100581599e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-3       3-1   -3.142728259861913e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-3       3-2   -4.140640161473054e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       3-3       3-3    8.815978581387224e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-0       3-0   -7.404741688720580e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-0       3-1    1.621239212646677e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-0       3-2   -1.170105244612809e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-0       3-3    3.962265087723471e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-0       3-0   -1.999272217942325e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-0       3-1    9.292339866625983e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-0       3-2   -2.768228497911302e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-0       3-3   -4.409689592288538e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-0       3-0    1.261566154344830e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-0       3-1    1.771078435193314e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-0       3-2    1.795899585901083e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-0       3-3   -3.144563100581598e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-0       3-0    5.541377852106689e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-0       3-1   -2.593702126290947e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-0       3-2    7.672631917264616e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-0       3-3    1.242623020384290e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-1       3-0    1.621239212646678e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-1       3-1   -2.253383717000577e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-1       3-2    1.607220219523158e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-1       3-3    1.995546275556730e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-1       3-0    9.292339866625636e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-1       3-1   -5.906173289163096e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-1       3-2    1.343175849379032e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-1       3-3    4.650006550235300e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-1       3-0    1.771078435193314e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-1       3-1    3.811566183504565e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-1       3-2    3.039539191280858e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-1       3-3   -3.142728259861845e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-1       3-0   -2.593702126290981e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-1       3-1    1.637004087231101e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-1       3-2   -3.751991732856030e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-1       3-3   -1.288817636549446e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-2       3-0   -1.170105244612809e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-2       3-1    1.607220219523158e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-2       3-2   -3.776084548214437e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-2       3-3   -1.181958275632975e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-2       3-0   -2.768228497911302e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-2       3-1    1.343175849379021e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-2       3-2   -9.600844695299749e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-2       3-3   -1.684755883826707e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-2       3-0    1.795899585901082e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-2       3-1    3.039539191280858e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-2       3-2    6.377867042855724e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-2       3-3   -4.140640161473054e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-2       3-0    7.672631917264616e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-2       3-1   -3.751991732855853e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-2       3-2    2.661030826383537e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-2       3-3    4.708358783568324e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-3       3-0    3.962265087723471e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-3       3-1    1.995546275556739e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-3       3-2   -1.181958275632975e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       3-3       3-3   -5.204268186958048e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-3       3-0   -4.409689592286318e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-3       3-1    4.650006550235299e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-3       3-2   -1.684755883826895e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       3-3       3-3   -1.293256617135360e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-3       3-0   -3.144563100581599e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-3       3-1   -3.142728259861870e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-3       3-2   -4.140640161473054e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       3-3       3-3    8.815978581387139e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-3       3-0    1.242623020386011e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-3       3-1   -1.288817636549446e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-3       3-2    4.708358783568400e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       3-3       3-3    3.584446870297007e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-0       4-0    3.470929748641858e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-0       4-1   -5.755658776922638e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-0       4-2    5.884383975180221e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-0       4-3    2.457437435059745e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-0       4-0    6.237734803067767e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-0       4-1    1.484377104577359e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-0       4-2    8.704039875421782e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-0       4-3   -6.766788868675382e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-0       4-0    4.834314628562037e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-0       4-1   -8.204806581694754e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-0       4-2    8.243990933619512e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-0       4-3    3.436251415748213e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-0       4-0    6.658538340747794e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-0       4-1   -2.574871607823913e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-0       4-2    9.289322479048287e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-0       4-3   -2.859619616477466e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-1       4-0   -5.755658776922638e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-1       4-1    9.638844595176898e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-1       4-2   -8.499753783882888e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-1       4-3   -8.997498329196587e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-1       4-0    1.484377104577359e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-1       4-1    1.856964364122009e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-1       4-2    3.848632540902145e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-1       4-3   -1.475004512213402e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-1       4-0   -8.204806581694754e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-1       4-1    1.342333602903778e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-1       4-2   -1.210607214414307e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-1       4-3   -1.260751176546798e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-1       4-0   -2.574871607823913e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-1       4-1    1.880984481568149e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-1       4-2   -3.078082530495410e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-1       4-3   -1.397363347372269e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-2       4-0    5.884383975180212e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-2       4-1   -8.499753783882888e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-2       4-2    1.457001409499597e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-2       4-3    1.079594004435839e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-2       4-0    8.704039875421782e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-2       4-1    3.848632540902145e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-2       4-2    3.041485009562201e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-2       4-3   -6.056026532693017e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-2       4-0    8.243990933619512e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-2       4-1   -1.210607214414307e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-2       4-2    2.029207791474160e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-2       4-3    1.536491486136132e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-2       4-0    9.289322479048287e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-2       4-1   -3.078082530495410e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-2       4-2    2.855827420063374e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-2       4-3    3.215456678678348e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-3       4-0    2.457437435059744e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-3       4-1   -8.997498329196587e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-3       4-2    1.079594004435839e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       4-3       4-3    1.894300923436015e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-3       4-0   -6.766788868675373e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-3       4-1   -1.475004512213402e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-3       4-2   -6.056026532693018e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       4-3       4-3    4.116248395369831e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-3       4-0    3.436251415748218e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-3       4-1   -1.260751176546796e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-3       4-2    1.536491486136132e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       4-3       4-3    2.640309114497233e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-3       4-0   -2.859619616477466e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-3       4-1   -1.397363347372269e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-3       4-2    3.215456678678348e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       4-3       4-3    3.592219261874056e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-0       4-0    6.237734803067881e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-0       4-1    1.484377104577440e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-0       4-2    8.704039875421805e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-0       4-3   -6.766788868675505e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-0       4-0    1.027291122418720e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-0       4-1   -1.730024151988892e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-0       4-2    1.748391502655766e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-0       4-3    7.292331303203656e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-0       4-0    8.686486339513980e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-0       4-1    2.127465200901815e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-0       4-2    1.212835427059640e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-0       4-3   -8.864620499426764e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-0       4-0   -8.165717827156087e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-0       4-1    1.428653828664232e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-0       4-2   -1.403456133657919e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-0       4-3   -5.834645287068095e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-1       4-0    1.484377104577440e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-1       4-1    1.856964364122024e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-1       4-2    3.848632540902286e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-1       4-3   -1.475004512213375e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-1       4-0   -1.730024151988892e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-1       4-1    2.852575638210560e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-1       4-2   -2.553360302988870e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-1       4-3   -2.673663113715486e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-1       4-0    2.127465200901815e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-1       4-1    2.587211931565825e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-1       4-2    5.349345699236160e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-1       4-3   -2.057652598521619e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-1       4-0    1.428653828664232e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-1       4-1   -2.266982260769828e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-1       4-2    2.105623473814782e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-1       4-3    2.146769740517087e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-2       4-0    8.704039875421805e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-2       4-1    3.848632540902285e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-2       4-2    3.041485009562181e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-2       4-3   -6.056026532693422e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-2       4-0    1.748391502655766e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-2       4-1   -2.553360302988870e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-2       4-2    4.312139055081928e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-2       4-3    3.241512883780941e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-2       4-0    1.212835427059640e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-2       4-1    5.349345699236160e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-2       4-2    4.240561938561446e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-2       4-3   -8.319418967461429e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-2       4-0   -1.403456133657919e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-2       4-1    2.105623473814783e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-2       4-2   -3.427334790444160e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-2       4-3   -2.669859780430970e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-3       4-0   -6.766788868675505e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-3       4-1   -1.475004512213375e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-3       4-2   -6.056026532693424e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       4-3       4-3    4.116248395369795e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-3       4-0    7.292331303203659e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-3       4-1   -2.673663113715487e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-3       4-2    3.241512883780940e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       4-3       4-3    5.609273665280424e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-3       4-0   -8.864620499426755e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-3       4-1   -2.057652598521619e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-3       4-2   -8.319418967461429e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       4-3       4-3    5.743186269856972e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-3       4-0   -5.834645287068098e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-3       4-1    2.146769740517087e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-3       4-2   -2.669859780430971e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       4-3       4-3   -4.464174378211299e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-0       4-0    4.834314628562039e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-0       4-1   -8.204806581694754e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-0       4-2    8.243990933619512e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-0       4-3    3.436251415748213e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-0       4-0    8.686486339514231e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-0       4-1    2.127465200901312e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-0       4-2    1.212835427059675e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-0       4-3   -8.864620499426976e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-0       4-0    1.674617460629549e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-0       4-1   -2.880745716788718e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-0       4-2    2.865616393028529e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-0       4-3    1.193059292430784e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-0       4-0   -1.044293288826940e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-0       4-1   -2.660984640943167e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-0       4-2   -1.458929453291039e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-0       4-3    1.004511817636994e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-1       4-0   -8.204806581694754e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-1       4-1    1.342333602903778e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-1       4-2   -1.210607214414307e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-1       4-3   -1.260751176546797e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-1       4-0    2.127465200901312e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-1       4-1    2.587211931565898e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-1       4-2    5.349345699235516e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-1       4-3   -2.057652598521669e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-1       4-0   -2.880745716788718e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-1       4-1    4.649536193904715e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-1       4-2   -4.248390493355161e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-1       4-3   -4.382802234815474e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-1       4-0   -2.660984640943167e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-1       4-1   -3.111939238700299e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-1       4-2   -6.472405836726411e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-1       4-3    2.478053037233067e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-2       4-0    8.243990933619512e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-2       4-1   -1.210607214414307e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-2       4-2    2.029207791474160e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-2       4-3    1.536491486136132e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-2       4-0    1.212835427059675e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-2       4-1    5.349345699235518e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-2       4-2    4.240561938561555e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-2       4-3   -8.319418967461327e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-2       4-0    2.865616393028529e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-2       4-1   -4.248390493355161e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-2       4-2    7.029014141349985e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-2       4-3    5.389690105931554e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-2       4-0   -1.458929453291039e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-2       4-1   -6.472405836726410e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-2       4-2   -5.104338938991206e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-2       4-3    9.940592757939942e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-3       4-0    3.436251415748219e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-3       4-1   -1.260751176546796e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-3       4-2    1.536491486136132e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       4-3       4-3    2.640309114497234e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-3       4-0   -8.864620499426974e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-3       4-1   -2.057652598521670e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-3       4-2   -8.319418967461328e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       4-3       4-3    5.743186269857087e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-3       4-0    1.193059292430785e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-3       4-1   -4.382802234815477e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-3       4-2    5.389690105931555e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       4-3       4-3    9.150052423011983e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-3       4-0    1.004511817636994e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-3       4-1    2.478053037233067e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-3       4-2    9.940592757939945e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       4-3       4-3   -6.918056107448498e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-0       4-0    6.658538340749663e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-0       4-1   -2.574871607824374e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-0       4-2    9.289322479059316e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-0       4-3   -2.859619616479586e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-0       4-0   -8.165717827156096e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-0       4-1    1.428653828664232e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-0       4-2   -1.403456133657919e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-0       4-3   -5.834645287068098e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-0       4-0   -1.044293288826883e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-0       4-1   -2.660984640944183e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-0       4-2   -1.458929453291038e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-0       4-3    1.004511817636268e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-0       4-0    2.264591343915777e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-0       4-1   -3.977084278520565e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-0       4-2    3.896038915580649e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-0       4-3    1.619112718130730e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-1       4-0   -2.574871607824374e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-1       4-1    1.880984481569874e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-1       4-2   -3.078082530495621e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-1       4-3   -1.397363347373844e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-1       4-0    1.428653828664232e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-1       4-1   -2.266982260769829e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-1       4-2    2.105623473814783e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-1       4-3    2.146769740517087e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-1       4-0   -2.660984640944182e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-1       4-1   -3.111939238700246e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-1       4-2   -6.472405836726598e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-1       4-3    2.478053037233157e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-1       4-0   -3.977084278520565e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-1       4-1    6.286874346605303e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-1       4-2   -5.860824096016946e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-1       4-3   -5.959669049730612e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-2       4-0    9.289322479059316e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-2       4-1   -3.078082530495621e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-2       4-2    2.855827420066223e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-2       4-3    3.215456678679340e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-2       4-0   -1.403456133657918e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-2       4-1    2.105623473814783e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-2       4-2   -3.427334790444160e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-2       4-3   -2.669859780430970e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-2       4-0   -1.458929453291038e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-2       4-1   -6.472405836726600e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-2       4-2   -5.104338938991262e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-2       4-3    9.940592757938635e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-2       4-0    3.896038915580648e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-2       4-1   -5.860824096016946e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-2       4-2    9.504931858269732e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-2       4-3    7.430434611799295e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-3       4-0   -2.859619616479586e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-3       4-1   -1.397363347373844e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-3       4-2    3.215456678679340e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       4-3       4-3    3.592219261876863e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-3       4-0   -5.834645287068098e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-3       4-1    2.146769740517087e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-3       4-2   -2.669859780430971e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       4-3       4-3   -4.464174378211299e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-3       4-0    1.004511817636268e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-3       4-1    2.478053037233157e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-3       4-2    9.940592757938637e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       4-3       4-3   -6.918056107448535e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-3       4-0    1.619112718130731e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-3       4-1   -5.959669049730618e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-3       4-2    7.430434611799294e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       4-3       4-3    1.238202952346929e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-0       5-0    1.617722841018002e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-0       5-1    1.611059558059760e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-0       5-2    2.193984799913066e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-0       5-3    6.745830104633744e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-0       5-0   -1.134455021605454e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-0       5-1    4.094138139913707e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-0       5-2   -1.612336244336864e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-0       5-3   -2.036290452676543e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-0       5-0    2.201496502678023e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-0       5-1    2.221490601923753e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-0       5-2    2.984819537362400e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-0       5-3    5.138425080252309e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-0       5-0    2.814365150745277e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-0       5-1   -9.515198156078695e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-0       5-2    3.887738110801330e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-0       5-3   -1.414505966591111e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-1       5-0    1.611059558059732e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-1       5-1    4.715388038434899e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-1       5-2    2.269532645620491e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-1       5-3   -3.583431430609223e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-1       5-0    4.094138139913707e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-1       5-1   -3.340553619094414e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-1       5-2    9.282316694136377e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-1       5-3    2.757680353027148e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-1       5-0    2.221490601923752e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-1       5-1    6.415731951627069e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-1       5-2    3.056976049775586e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-1       5-3   -4.873053314266748e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-1       5-0   -9.515198156078695e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-1       5-1    8.144637745711195e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-1       5-2   -1.100250442897302e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-1       5-3   -6.132527105838997e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-2       5-0    2.193984799913066e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-2       5-1    2.269532645620491e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-2       5-2    7.520893114168553e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-2       5-3   -2.756298518435941e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-2       5-0   -1.612336244336313e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-2       5-1    9.282316694136377e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-2       5-2   -5.472503048264271e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-2       5-3   -1.517884088121525e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-2       5-0    2.984819537362400e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-2       5-1    3.056976049775364e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-2       5-2    1.022999737497274e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-2       5-3   -3.632289712177758e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-2       5-0    3.887738110797172e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-2       5-1   -1.100250442897302e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-2       5-2    1.276536406274657e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-2       5-3    1.068454938351915e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-3       5-0    6.745830104688821e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-3       5-1   -3.583431430609223e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-3       5-2   -2.756298518435722e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-0       5-3       5-3    9.939166311457869e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-3       5-0   -2.036290452676543e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-3       5-1    2.757680353031519e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-3       5-2   -1.517884088121525e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-1       5-3       5-3   -7.542259557922204e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-3       5-0    5.138425080274513e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-3       5-1   -4.873053314266748e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-3       5-2   -3.632289712178639e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-2       5-3       5-3    1.351546530486151e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-3       5-0   -1.414505966591112e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-3       5-1   -6.132527105839039e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-3       5-2    1.068454938351915e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-0       2-3       5-3       5-3    1.655218105604418e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-0       5-0   -1.134455021606009e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-0       5-1    4.094138139913706e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-0       5-2   -1.612336244336864e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-0       5-3   -2.036290452676543e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-0       5-0    4.715207194820414e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-0       5-1    4.843491491138168e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-0       5-2    6.393595888610783e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-0       5-3    2.626639156041999e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-0       5-0   -1.004531796611225e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-0       5-1    7.386442556116529e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-0       5-2   -1.478728091014476e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-0       5-3   -2.623396695851560e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-0       5-0   -3.604236994532847e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-0       5-1   -3.825242305600846e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-0       5-2   -4.884599643762295e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-0       5-3   -8.715839631308515e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-1       5-0    4.094138139913706e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-1       5-1   -3.340553619098855e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-1       5-2    9.282316694136377e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-1       5-3    2.757680353027148e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-1       5-0    4.843491491137307e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-1       5-1    1.374225930995464e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-1       5-2    6.812241860955625e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-1       5-3   -1.043976565488681e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-1       5-0    7.386442556116528e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-1       5-1   -2.966227648633600e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-1       5-2    1.494045481578634e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-1       5-3    2.553969372563053e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-1       5-0   -3.825242305600957e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-1       5-1   -1.050074880820430e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-1       5-2   -5.265413616769970e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-1       5-3    7.969886560461473e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-2       5-0   -1.612336244336313e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-2       5-1    9.282316694136377e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-2       5-2   -5.472503048255389e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-2       5-3   -1.517884088121525e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-2       5-0    6.393595888610784e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-2       5-1    6.812241860957843e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-2       5-2    2.191439406717343e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-2       5-3   -8.261510033917063e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-2       5-0   -1.478728091016141e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-2       5-1    1.494045481578634e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-2       5-2   -4.897300549957656e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-2       5-3   -2.318264786689687e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-2       5-0   -4.884599643762295e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-2       5-1   -5.265413616769956e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-2       5-2   -1.673689264349222e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-2       5-3    6.258192856139363e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-3       5-0   -2.036290452676543e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-3       5-1    2.757680353031519e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-3       5-2   -1.517884088121525e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-0       5-3       5-3   -7.542259557922204e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-3       5-0    2.626639155975993e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-3       5-1   -1.043976565488681e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-3       5-2   -8.261510033917299e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-1       5-3       5-3    2.895528146981674e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-3       5-0   -2.623396695851560e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-3       5-1    2.553969372562429e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-3       5-2   -2.318264786689687e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-2       5-3       5-3   -6.856067558244083e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-3       5-0   -8.715839631342516e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-3       5-1    7.969886560461473e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-3       5-2    6.258192856137809e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-1       2-3       5-3       5-3   -2.210305087626455e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-0       5-0    2.201496502678023e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-0       5-1    2.221490601923753e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-0       5-2    2.984819537362400e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-0       5-3    5.138425080274513e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-0       5-0   -1.004531796615701e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-0       5-1    7.386442556116529e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-0       5-2   -1.478728091014753e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-0       5-3   -2.623396695851560e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-0       5-0    7.523838056089197e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-0       5-1    8.026461478826831e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-0       5-2    1.019904987748915e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-0       5-3    5.446556527117195e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-0       5-0    7.722582614600615e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-0       5-1   -1.071092013374094e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-0       5-2    1.202437146129295e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-0       5-3    2.915386308686165e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-1       5-0    2.221490601923752e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-1       5-1    6.415731951627067e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-1       5-2    3.056976049775586e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-1       5-3   -4.873053314266747e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-1       5-0    7.386442556116528e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-1       5-1   -2.966227648643037e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-1       5-2    1.494045481578634e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-1       5-3    2.553969372565273e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-1       5-0    8.026461478826611e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-1       5-1    2.192376682253524e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-1       5-2    1.126903814655294e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-1       5-3   -1.664678676743983e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-1       5-0   -1.071092013374094e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-1       5-1    2.287784771198651e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-1       5-2   -2.013653178739073e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-1       5-3   -2.093730181448733e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-2       5-0    2.984819537362400e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-2       5-1    3.056976049775364e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-2       5-2    1.022999737497274e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-2       5-3   -3.632289712178646e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-2       5-0   -1.478728091016696e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-2       5-1    1.494045481578634e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-2       5-2   -4.897300549962097e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-2       5-3   -2.318264786689687e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-2       5-0    1.019904987748915e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-2       5-1    1.126903814655249e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-2       5-2    3.495174901063422e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-2       5-3   -1.364443954306446e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-2       5-0    1.202437146130128e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-2       5-1   -2.013653178739074e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-2       5-2    3.811484729166620e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-2       5-3    3.006055304578511e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-3       5-0    5.138425080274513e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-3       5-1   -4.873053314266747e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-3       5-2   -3.632289712178639e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-0       5-3       5-3    1.351546530486151e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-3       5-0   -2.623396695851560e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-3       5-1    2.553969372564649e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-3       5-2   -2.318264786689687e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-1       5-3       5-3   -6.856067558252965e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-3       5-0    5.446556527117195e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-3       5-1   -1.664678676743983e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-3       5-2   -1.364443954306535e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-2       5-3       5-3    4.616867654011791e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-3       5-0    2.915386308686165e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-3       5-1   -2.093730181475933e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-3       5-2    3.006055304578511e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-2       2-3       5-3       5-3    5.446978171046446e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-0       5-0    2.814365150736823e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-0       5-1   -9.515198156078692e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-0       5-2    3.887738110792047e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-0       5-3   -1.414505966591111e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-0       5-0   -3.604236994532847e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-0       5-1   -3.825242305600846e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-0       5-2   -4.884599643762295e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-0       5-3   -8.715839631352924e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-0       5-0    7.722582614645023e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-0       5-1   -1.071092013374094e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-0       5-2    1.202437146129295e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-0       5-3    2.915386308686165e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-0       5-0    9.962632967736461e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-0       5-1    1.095597330440650e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-0       5-2    1.350097838891158e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-0       5-3    9.556737040574553e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-1       5-0   -9.515198156078692e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-1       5-1    8.144637745695203e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-1       5-2   -1.100250442897302e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-1       5-3   -6.132527105819908e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-1       5-0   -3.825242305600957e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-1       5-1   -1.050074880820430e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-1       5-2   -5.265413616769970e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-1       5-3    7.969886560461474e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-1       5-0   -1.071092013374094e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-1       5-1    2.287784771216415e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-1       5-2   -2.013653178739073e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-1       5-3   -2.093730181466497e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-1       5-0    1.095597330440600e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-1       5-1    2.902451273746276e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-1       5-2    1.534430726066491e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-1       5-3   -2.202693364923135e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-2       5-0    3.887738110790703e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-2       5-1   -1.100250442897302e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-2       5-2    1.276536406268664e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-2       5-3    1.068454938351915e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-2       5-0   -4.884599643762295e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-2       5-1   -5.265413616770858e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-2       5-2   -1.673689264349222e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-2       5-3    6.258192856137643e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-2       5-0    1.202437146130128e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-2       5-1   -2.013653178739073e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-2       5-2    3.811484729166620e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-2       5-3    3.006055304578511e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-2       5-0    1.350097838891158e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-2       5-1    1.534430726066491e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-2       5-2    4.625901504443000e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-2       5-3   -1.853685192236101e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-3       5-0   -1.414505966591111e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-3       5-1   -6.132527105817134e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-3       5-2    1.068454938351915e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-0       5-3       5-3    1.655218105600160e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-3       5-0   -8.715839631341822e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-3       5-1    7.969886560461474e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-3       5-2    6.258192856137865e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-1       5-3       5-3   -2.210305087626455e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-3       5-0    2.915386308686165e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-3       5-1   -2.093730181475933e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-3       5-2    3.006055304578511e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-2       5-3       5-3    5.446978171081973e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-3       5-0    9.556737040580104e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-3       5-1   -2.202693364923135e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-3       5-2   -1.853685192236057e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "2-3       2-3       5-3       5-3    6.108722081447472e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-0       4-0    1.073138700123088e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-0       4-1   -1.488016834918449e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-0       4-2    1.668144438558419e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-0       4-3   -4.123484389272300e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-0       4-0    1.127059324486421e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-0       4-1    1.971772317771408e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-0       4-2    1.510954530939801e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-0       4-3    1.403456635632403e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-0       4-0    1.292814027575350e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-0       4-1   -2.070080904036721e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-0       4-2    2.055758128192022e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-0       4-3   -5.722627492972250e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-0       4-0   -2.627699874681039e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-0       4-1   -2.790315940318596e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-0       4-2   -2.980309819085811e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-0       4-3   -1.532312671526572e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-1       4-0   -1.488016834918449e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-1       4-1    2.875703475986870e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-1       4-2   -2.016726792188500e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-1       4-3   -2.288492812212910e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-1       4-0    1.971772317771408e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-1       4-1    3.265237981211911e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-1       4-2    1.406634271613977e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-1       4-3   -2.453494277607242e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-1       4-0   -2.070080904036721e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-1       4-1    3.426865830520982e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-1       4-2   -2.807645537110748e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-1       4-3   -2.765352048970703e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-1       4-0   -2.790315940318596e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-1       4-1   -7.861688059796803e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-1       4-2   -3.692835667047241e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-1       4-3    5.272847076514029e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-2       4-0    1.668144438558419e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-2       4-1   -2.016726792188500e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-2       4-2    4.043386833746361e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-2       4-3    2.369412567182496e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-2       4-0    1.510954530939800e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-2       4-1    1.406634271613977e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-2       4-2    5.178809385038964e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-2       4-3   -1.305896201418611e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-2       4-0    2.055758128192021e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-2       4-1   -2.807645537110748e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-2       4-2    4.730258657771981e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-2       4-3    3.300468769206730e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-2       4-0   -2.980309819085810e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-2       4-1   -3.692835667047241e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-2       4-2   -1.294060033899465e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-2       4-3    4.148863446501886e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-3       4-0   -4.123484389272293e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-3       4-1   -2.288492812212909e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-3       4-2    2.369412567182496e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       4-3       4-3    4.757807284664856e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-3       4-0    1.403456635632403e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-3       4-1   -2.453494277607242e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-3       4-2   -1.305896201418613e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       4-3       4-3    6.821444079312564e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-3       4-0   -5.722627492972246e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-3       4-1   -2.765352048970703e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-3       4-2    3.300468769206729e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       4-3       4-3    5.460404696551785e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-3       4-0   -1.532312671526572e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-3       4-1    5.272847076514029e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-3       4-2    4.148863446501886e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       4-3       4-3   -1.744258979665496e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-0       4-0    1.127059324486452e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-0       4-1    1.971772317771358e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-0       4-2    1.510954530939787e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-0       4-3    1.403456635632227e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-0       4-0    2.897454764693502e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-0       4-1   -4.409611532088001e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-0       4-2    4.569142571786174e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-0       4-3   -1.220005087537026e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-0       4-0    1.656170445459249e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-0       4-1    4.275610908153691e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-0       4-2    2.216123504841645e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-0       4-3    2.152660983131409e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-0       4-0   -1.746238860004792e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-0       4-1    3.498590326914314e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-0       4-2   -2.892983263007729e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-0       4-3    9.637278589473506e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-1       4-0    1.971772317771358e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-1       4-1    3.265237981211924e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-1       4-2    1.406634271614202e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-1       4-3   -2.453494277607149e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-1       4-0   -4.409611532088001e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-1       4-1    7.711360078936286e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-1       4-2   -5.979281384312491e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-1       4-3   -6.190749040250628e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-1       4-0    4.275610908153691e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-1       4-1    4.797351595874407e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-1       4-2    3.852333248055048e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-1       4-3   -3.596886430309106e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-1       4-0    3.498590326914314e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-1       4-1   -4.532962028896254e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-1       4-2    4.749676936318994e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-1       4-3    3.754323428364159e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-2       4-0    1.510954530939787e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-2       4-1    1.406634271614202e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-2       4-2    5.178809385038873e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-2       4-3   -1.305896201419036e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-2       4-0    4.569142571786179e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-2       4-1   -5.979281384312491e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-2       4-2    1.071811346941222e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-2       4-3    7.027513147085415e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-2       4-0    2.216123504841645e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-2       4-1    3.852333248055045e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-2       4-2    7.605702512938781e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-2       4-3   -3.707711328508712e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-2       4-0   -2.892983263007729e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-2       4-1    4.749676936318993e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-2       4-2   -6.029669131463461e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-2       4-3   -5.587456981729408e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-3       4-0    1.403456635632228e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-3       4-1   -2.453494277607149e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-3       4-2   -1.305896201419040e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       4-3       4-3    6.821444079312407e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-3       4-0   -1.220005087537023e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-3       4-1   -6.190749040250626e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-3       4-2    7.027513147085415e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       4-3       4-3    1.246263785417542e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-3       4-0    2.152660983131409e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-3       4-1   -3.596886430309106e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-3       4-2   -3.707711328508713e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       4-3       4-3    1.001269610338135e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-3       4-0    9.637278589473506e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-3       4-1    3.754323428364160e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-3       4-2   -5.587456981729409e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       4-3       4-3   -6.682462493709195e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-0       4-0    1.292814027575352e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-0       4-1   -2.070080904036721e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-0       4-2    2.055758128192021e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-0       4-3   -5.722627492972255e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-0       4-0    1.656170445459168e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-0       4-1    4.275610908153883e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-0       4-2    2.216123504841552e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-0       4-3    2.152660983131611e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-0       4-0    4.103548479739345e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-0       4-1   -7.198511527673288e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-0       4-2    6.628890539755143e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-0       4-3   -1.986780466139075e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-0       4-0   -2.077898651538590e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-0       4-1   -6.437915935010555e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-0       4-2   -2.778429422029726e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-0       4-3   -2.739757211221082e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-1       4-0   -2.070080904036721e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-1       4-1    3.426865830520986e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-1       4-2   -2.807645537110748e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-1       4-3   -2.765352048970703e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-1       4-0    4.275610908153880e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-1       4-1    4.797351595874198e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-1       4-2    3.852333248054902e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-1       4-3   -3.596886430309050e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-1       4-0   -7.198511527673290e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-1       4-1    1.079138745180918e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-1       4-2   -9.767412326336121e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-1       4-3   -8.793896213649845e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-1       4-0   -6.437915935010554e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-1       4-1   -6.020084189686577e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-1       4-2   -6.278719075328843e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-1       4-3    4.510441213415232e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-2       4-0    2.055758128192021e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-2       4-1   -2.807645537110748e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-2       4-2    4.730258657771982e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-2       4-3    3.300468769206730e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-2       4-0    2.216123504841552e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-2       4-1    3.852333248054900e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-2       4-2    7.605702512938568e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-2       4-3   -3.707711328508291e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-2       4-0    6.628890539755143e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-2       4-1   -9.767412326336121e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-2       4-2    1.469198193627501e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-2       4-3    1.148554554598138e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-2       4-0   -2.778429422029726e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-2       4-1   -6.278719075328845e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-2       4-2   -9.545144792201227e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-2       4-3    6.120708370116560e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-3       4-0   -5.722627492972251e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-3       4-1   -2.765352048970703e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-3       4-2    3.300468769206729e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       4-3       4-3    5.460404696551787e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-3       4-0    2.152660983131612e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-3       4-1   -3.596886430309050e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-3       4-2   -3.707711328508292e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       4-3       4-3    1.001269610338122e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-3       4-0   -1.986780466139075e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-3       4-1   -8.793896213649838e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-3       4-2    1.148554554598138e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       4-3       4-3    1.671055473304539e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-3       4-0   -2.739757211221082e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-3       4-1    4.510441213415232e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-3       4-2    6.120708370116557e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       4-3       4-3   -1.256558058591073e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-0       4-0   -2.627699874680513e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-0       4-1   -2.790315940318496e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-0       4-2   -2.980309819085187e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-0       4-3   -1.532312671527525e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-0       4-0   -1.746238860004789e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-0       4-1    3.498590326914314e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-0       4-2   -2.892983263007729e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-0       4-3    9.637278589473504e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-0       4-0   -2.077898651538760e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-0       4-1   -6.437915935008838e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-0       4-2   -2.778429422029861e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-0       4-3   -2.739757211221062e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-0       4-0    4.746792230191796e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-0       4-1   -9.735420397280218e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-0       4-2    7.899823788741300e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-0       4-3   -2.679876513827722e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-1       4-0   -2.790315940318496e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-1       4-1   -7.861688059795357e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-1       4-2   -3.692835667046932e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-1       4-3    5.272847076515497e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-1       4-0    3.498590326914314e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-1       4-1   -4.532962028896249e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-1       4-2    4.749676936318993e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-1       4-3    3.754323428364161e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-1       4-0   -6.437915935008837e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-1       4-1   -6.020084189686931e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-1       4-2   -6.278719075326553e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-1       4-3    4.510441213415450e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-1       4-0   -9.735420397280218e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-1       4-1    1.228921193001857e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-1       4-2   -1.321812440466688e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-1       4-3   -1.020633135937038e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-2       4-0   -2.980309819085187e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-2       4-1   -3.692835667046932e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-2       4-2   -1.294060033899537e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-2       4-3    4.148863446501435e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-2       4-0   -2.892983263007729e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-2       4-1    4.749676936318993e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-2       4-2   -6.029669131463459e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-2       4-3   -5.587456981729408e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-2       4-0   -2.778429422029862e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-2       4-1   -6.278719075326556e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-2       4-2   -9.545144792201744e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-2       4-3    6.120708370117314e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-2       4-0    7.899823788741313e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-2       4-1   -1.321812440466688e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-2       4-2    1.626781491457168e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-2       4-3    1.555087070078268e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-3       4-0   -1.532312671527525e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-3       4-1    5.272847076515497e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-3       4-2    4.148863446501435e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       4-3       4-3   -1.744258979665776e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-3       4-0    9.637278589473484e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-3       4-1    3.754323428364161e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-3       4-2   -5.587456981729409e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       4-3       4-3   -6.682462493709195e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-3       4-0   -2.739757211221062e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-3       4-1    4.510441213415450e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-3       4-2    6.120708370117313e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       4-3       4-3   -1.256558058591125e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-3       4-0   -2.679876513827722e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-3       4-1   -1.020633135937038e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-3       4-2    1.555087070078269e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       4-3       4-3    1.792830579721161e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-0       5-0    1.049966686872319e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-0       5-1    2.118320086854757e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-0       5-2    1.398149841246068e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-0       5-3    1.794732346634240e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-0       5-0   -2.782318077043755e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-0       5-1    1.626692440388032e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-0       5-2   -4.177681715400278e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-0       5-3    2.162390724183759e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-0       5-0    1.455063411919339e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-0       5-1    3.160954223027196e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-0       5-2    1.937205934971447e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-0       5-3    2.331790607338190e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-0       5-0   -8.382022392351173e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-0       5-1   -1.491048406886653e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-0       5-2   -1.213806298468087e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-0       5-3    4.949516436757437e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-1       5-0    2.118320086854840e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-1       5-1    3.023959028592518e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-1       5-2    2.686390845497651e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-1       5-3   -2.223531418966234e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-1       5-0    1.626692440388031e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-1       5-1   -8.654210166057992e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-1       5-2    1.925227682122984e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-1       5-3    7.701903028187161e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-1       5-0    3.160954223027189e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-1       5-1    4.190127062484016e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-1       5-2    4.067731979551479e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-1       5-3   -3.079913884294596e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-1       5-0   -1.491048406886652e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-1       5-1   -2.577163971067127e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-1       5-2   -2.113802377974285e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-1       5-3    2.209789360219151e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-2       5-0    1.398149841246068e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-2       5-1    2.686390845497651e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-2       5-2    4.738437721261922e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-2       5-3   -2.928932744898347e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-2       5-0   -4.177681715399990e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-2       5-1    1.925227682122984e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-2       5-2   -1.505439561271757e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-2       5-3   -2.024898048489133e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-2       5-0    1.937205934971447e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-2       5-1    4.067731979551534e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-2       5-2    6.564523881102656e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-2       5-3   -4.503545530141519e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-2       5-0   -1.213806298468104e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-2       5-1   -2.113802377974286e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-2       5-2   -4.424113085599863e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-2       5-3    2.586832517557065e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-3       5-0    1.794732346634049e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-3       5-1   -2.223531418966234e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-3       5-2   -2.928932744898513e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-0       5-3       5-3    6.146622414961761e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-3       5-0    2.162390724183759e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-3       5-1    7.701903028185825e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-3       5-2   -2.024898048489133e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-1       5-3       5-3   -2.161294288762829e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-3       5-0    2.331790607337757e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-3       5-1   -3.079913884294596e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-3       5-2   -4.503545530141879e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-2       5-3       5-3    8.513625150935836e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-3       5-0    4.949516436757437e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-3       5-1    2.209789360219155e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-3       5-2    2.586832517557065e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-0       3-3       5-3       5-3   -6.260704733348534e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-0       5-0   -2.782318077043337e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-0       5-1    1.626692440388031e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-0       5-2   -4.177681715400553e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-0       5-3    2.162390724183760e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-0       5-0    3.103551973656657e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-0       5-1    6.722999944908763e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-0       5-2    4.132200320632317e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-0       5-3    5.271487431445770e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-0       5-0   -1.684694662965972e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-0       5-1    4.857799912321841e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-0       5-2   -7.243591252869598e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-0       5-3    2.966876000380496e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-0       5-0   -2.446353695287843e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-0       5-1   -5.765914994335875e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-0       5-2   -3.256146881492778e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-0       5-3   -3.882918718372142e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-1       5-0    1.626692440388031e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-1       5-1   -8.654210166061336e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-1       5-2    1.925227682122984e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-1       5-3    7.701903028186675e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-1       5-0    6.722999944908770e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-1       5-1    8.937640936347176e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-1       5-2    8.598028267082012e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-1       5-3   -6.570326654781184e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-1       5-0    4.857799912321841e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-1       5-1   -1.110526559194912e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-1       5-2    6.349168454230393e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-1       5-3    2.166563210893166e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-1       5-0   -5.765914994335764e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-1       5-1   -7.043576008660644e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-1       5-2   -7.482886273359747e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-1       5-3    5.174907470441600e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-2       5-0   -4.177681715400822e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-2       5-1    1.925227682122984e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-2       5-2   -1.505439561271878e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-2       5-3   -2.024898048489133e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-2       5-0    4.132200320632317e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-2       5-1    8.598028267082345e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-2       5-2    1.400319221540325e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-2       5-3   -9.465374431224749e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-2       5-0   -7.243591252870975e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-2       5-1    6.349168454230393e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-2       5-2   -3.187353646077109e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-2       5-3   -7.302608516173760e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-2       5-0   -3.256146881492777e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-2       5-1   -7.482886273359720e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-2       5-2   -1.103218355769356e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-2       5-3    8.363834304514560e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-3       5-0    2.162390724183760e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-3       5-1    7.701903028186953e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-3       5-2   -2.024898048489133e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-0       5-3       5-3   -2.161294288762621e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-3       5-0    5.271487431444104e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-3       5-1   -6.570326654781184e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-3       5-2   -9.465374431224749e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-1       5-3       5-3    1.816221282971263e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-3       5-0    2.966876000380496e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-3       5-1    2.166563210893096e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-3       5-2   -7.302608516173760e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-2       5-3       5-3   -6.183135287751579e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-3       5-0   -3.882918718375334e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-3       5-1    5.174907470441600e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-3       5-2    8.363834304514338e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-1       3-3       5-3       5-3   -1.430395454285890e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-0       5-0    1.455063411919339e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-0       5-1    3.160954223026974e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-0       5-2    1.937205934971447e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-0       5-3    2.331790607337080e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-0       5-0   -1.684694662968791e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-0       5-1    4.857799912321842e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-0       5-2   -7.243591252870292e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-0       5-3    2.966876000380496e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-0       5-0    5.048229334845544e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-0       5-1    1.175903812362495e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-0       5-2    6.720257169869505e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-0       5-3    8.490090319087963e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-0       5-0   -3.447120018705575e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-0       5-1   -8.662317564850452e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-0       5-2   -4.164790198874062e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-0       5-3   -3.523081025144189e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-1       5-0    3.160954223026745e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-1       5-1    4.190127062484016e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-1       5-2    4.067731979551479e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-1       5-3   -3.079913884294596e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-1       5-0    4.857799912321842e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-1       5-1   -1.110526559198243e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-1       5-2    6.349168454230393e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-1       5-3    2.166563210897884e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-1       5-0    1.175903812362415e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-1       5-1    1.453629069161441e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-1       5-2    1.516268972629038e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-1       5-3   -1.068262279686127e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-1       5-0   -8.662317564850452e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-1       5-1   -9.471913055710667e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-1       5-2   -1.162540566734104e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-1       5-3    5.911065993320173e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-2       5-0    1.937205934971447e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-2       5-1    4.067731979551534e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-2       5-2    6.564523881102656e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-2       5-3   -4.503545530141519e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-2       5-0   -7.243591252904281e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-2       5-1    6.349168454230393e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-2       5-2   -3.187353646077109e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-2       5-3   -7.302608516173760e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-2       5-0    6.720257169869505e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-2       5-1    1.516268972629061e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-2       5-2    2.277105620441823e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-2       5-3   -1.684661230076050e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-2       5-0   -4.164790198871858e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-2       5-1   -1.162540566734104e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-2       5-2   -1.380229784797564e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-2       5-3    1.365769501373846e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-3       5-0    2.331790607337757e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-3       5-1   -3.079913884294596e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-3       5-2   -4.503545530141879e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-0       5-3       5-3    8.513625150935837e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-3       5-0    2.966876000380496e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-3       5-1    2.166563210897815e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-3       5-2   -7.302608516173759e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-1       5-3       5-3   -6.183135287770453e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-3       5-0    8.490090319089351e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-3       5-1   -1.068262279686127e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-3       5-2   -1.684661230075674e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-2       5-3       5-3    2.952869017098397e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-3       5-0   -3.523081025144190e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-3       5-1    5.911065993317831e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-3       5-2    1.365769501373846e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-2       3-3       5-3       5-3   -1.638592067444905e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-0       5-0   -8.382022392351168e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-0       5-1   -1.491048406886653e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-0       5-2   -1.213806298468082e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-0       5-3    4.949516436757435e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-0       5-0   -2.446353695287843e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-0       5-1   -5.765914994335876e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-0       5-2   -3.256146881492777e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-0       5-3   -3.882918718375438e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-0       5-0   -3.447120018705582e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-0       5-1   -8.662317564850452e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-0       5-2   -4.164790198865301e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-0       5-3   -3.523081025144189e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-0       5-0    6.802878466452983e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-0       5-1    1.659760177473565e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-0       5-2    9.054581762865338e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-0       5-3    1.128535503683690e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-1       5-0   -1.491048406886652e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-1       5-1   -2.577163971064901e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-1       5-2   -2.113802377974285e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-1       5-3    2.209789360219107e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-1       5-0   -5.765914994335543e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-1       5-1   -7.043576008660643e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-1       5-2   -7.482886273359754e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-1       5-3    5.174907470441600e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-1       5-0   -8.662317564850452e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-1       5-1   -9.471913055710632e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-1       5-2   -1.162540566734104e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-1       5-3    5.911065993319756e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-1       5-0    1.659760177473575e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-1       5-1    1.958667656199152e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-1       5-2    2.152058150695352e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-1       5-3   -1.438975392759437e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-2       5-0   -1.213806298468099e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-2       5-1   -2.113802377974286e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-2       5-2   -4.424113085599824e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-2       5-3    2.586832517557065e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-2       5-0   -3.256146881492777e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-2       5-1   -7.482886273359720e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-2       5-2   -1.103218355769356e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-2       5-3    8.363834304514560e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-2       5-0   -4.164790198867522e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-2       5-1   -1.162540566734104e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-2       5-2   -1.380229784797647e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-2       5-3    1.365769501373846e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-2       5-0    9.054581762865338e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-2       5-1    2.152058150695297e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-2       5-2    3.067750820748340e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-2       5-3   -2.405443770409588e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-3       5-0    4.949516436757435e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-3       5-1    2.209789360219112e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-3       5-2    2.586832517557065e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-0       5-3       5-3   -6.260704733348325e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-3       5-0   -3.882918718376410e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-3       5-1    5.174907470441600e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-3       5-2    8.363834304514393e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-1       5-3       5-3   -1.430395454285890e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-3       5-0   -3.523081025144190e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-3       5-1    5.911065993318108e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-3       5-2    1.365769501373846e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-2       5-3       5-3   -1.638592067445016e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-3       5-0    1.128535503683857e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-3       5-1   -1.438975392759437e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-3       5-2   -2.405443770409765e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "3-3       3-3       5-3       5-3    3.977452917701449e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-0       5-0   -2.749812630600969e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-0       5-1   -3.078043062858376e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-0       5-2   -4.084199733828908e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-0       5-3   -2.313907770799054e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-0       5-0   -3.155086491351444e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-0       5-1   -2.374856592168591e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-0       5-2   -4.337833880511357e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-0       5-3   -2.822793054407524e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-0       5-0   -3.939719595486457e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-0       5-1   -4.359311345872023e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-0       5-2   -5.850895297076978e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-0       5-3   -3.163029162804693e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-0       5-0    7.521537027004103e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-0       5-1    1.220827707549312e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-0       5-2    1.054749127222759e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-0       5-3    1.091629478004341e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-1       5-0   -3.078043062858434e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-1       5-1   -8.516805712869846e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-1       5-2   -3.949325756563394e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-1       5-3    7.479021400317640e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-1       5-0   -2.374856592168509e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-1       5-1   -9.279727433625224e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-1       5-2   -2.880571058585414e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-1       5-3    7.224009491766846e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-1       5-0   -4.359311345872023e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-1       5-1   -1.220132850613081e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-1       5-2   -5.612905371621611e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-1       5-3    1.071094448104043e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-1       5-0    1.220827707549319e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-1       5-1    2.241395247038476e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-1       5-2    1.536318779421031e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-1       5-3   -1.799547610888852e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-2       5-0   -4.084199733828908e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-2       5-1   -3.949325756563371e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-2       5-2   -1.472666131343781e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-2       5-3    4.365606431797130e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-2       5-0   -4.337833880511357e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-2       5-1   -2.880571058585359e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-2       5-2   -1.499647882205514e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-2       5-3    3.024954473553457e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-2       5-0   -5.850895297076978e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-2       5-1   -5.612905371621604e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-2       5-2   -2.109337056927749e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-2       5-3    6.227068117149322e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-2       5-0    1.054749127222759e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-2       5-1    1.536318779421028e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-2       5-2    3.683955690278245e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-2       5-3   -1.657379041409152e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-3       5-0   -2.313907770798880e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-3       5-1    7.479021400317641e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-3       5-2    4.365606431797158e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-0       5-3       5-3   -2.098317228636979e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-3       5-0   -2.822793054406900e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-3       5-1    7.224009491766844e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-3       5-2    3.024954473552791e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-1       5-3       5-3   -2.009001094854225e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-3       5-0   -3.163029162804382e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-3       5-1    1.071094448104043e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-3       5-2    6.227068117149306e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-2       5-3       5-3   -3.004508733275824e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-3       5-0    1.091629478004338e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-3       5-1   -1.799547610888852e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-3       5-2   -1.657379041409158e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-0       4-3       5-3       5-3    5.008884391539085e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-0       5-0   -3.155086491351444e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-0       5-1   -2.374856592168813e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-0       5-2   -4.337833880511357e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-0       5-3   -2.822793054407524e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-0       5-0   -9.475285820729200e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-0       5-1   -9.315579691927348e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-0       5-2   -1.394342298019837e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-0       5-3   -6.880517675596832e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-0       5-0   -4.606555295141017e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-0       5-1   -5.740963848105020e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-0       5-2   -6.339931561521386e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-0       5-3   -6.070004790510912e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-0       5-0    8.693532935491385e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-0       5-1    7.521484092497883e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-0       5-2    1.271104334664337e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-0       5-3    5.286945153288017e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-1       5-0   -2.374856592168731e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-1       5-1   -9.279727433625224e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-1       5-2   -2.880571058586303e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-1       5-3    7.224009491766846e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-1       5-0   -9.315579691927348e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-1       5-1   -2.916356537663437e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-1       5-2   -1.197399942490950e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-1       5-3    2.526224971414817e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-1       5-0   -5.740963848105048e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-1       5-1   -1.355799190977935e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-1       5-2   -7.080820461285156e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-1       5-3    1.057221006713734e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-1       5-0    7.521484092497820e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-1       5-1    2.664155641631492e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-1       5-2    9.714497984921451e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-1       5-3   -2.285369644867212e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-2       5-0   -4.337833880511357e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-2       5-1   -2.880571058586247e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-2       5-2   -1.499647882205514e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-2       5-3    3.024954473553457e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-2       5-0   -1.394342298019837e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-2       5-1   -1.197399942490950e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-2       5-2   -5.003286746526438e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-2       5-3    1.326199441242662e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-2       5-0   -6.339931561521386e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-2       5-1   -7.080820461285142e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-2       5-2   -2.193037679558781e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-2       5-3    7.533896335955859e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-2       5-0    1.271104334664337e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-2       5-1    9.714497984921558e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-2       5-2    4.545162992418934e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-2       5-3   -1.081392657106269e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-3       5-0   -2.822793054409120e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-3       5-1    7.224009491766844e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-3       5-2    3.024954473552791e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-0       5-3       5-3   -2.009001094854225e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-3       5-0   -6.880517675597370e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-3       5-1    2.526224971414817e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-3       5-2    1.326199441242707e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-1       5-3       5-3   -7.080217118011683e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-3       5-0   -6.070004790505153e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-3       5-1    1.057221006713734e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-3       5-2    7.533896335956539e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-2       5-3       5-3   -2.940388320173311e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-3       5-0    5.286945153288019e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-3       5-1   -2.285369644867212e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-3       5-2   -1.081392657106271e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-1       4-3       5-3       5-3    6.399605663311412e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-0       5-0   -3.939719595486458e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-0       5-1   -4.359311345872066e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-0       5-2   -5.850895297076977e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-0       5-3   -3.163029162804832e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-0       5-0   -4.606555295141017e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-0       5-1   -5.740963848105048e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-0       5-2   -6.339931561521386e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-0       5-3   -6.070004790518684e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-0       5-0   -1.824860505578433e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-0       5-1   -1.549659949075027e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-0       5-2   -2.662147958427581e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-0       5-3   -1.119970267389207e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-0       5-0    5.777752975611131e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-0       5-1    8.864520385006191e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-0       5-2    7.958531960397842e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-0       5-3    9.399142261971079e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-1       5-0   -4.359311345872073e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-1       5-1   -1.220132850613081e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-1       5-2   -5.612905371621652e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-1       5-3    1.071094448104043e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-1       5-0   -5.740963848105076e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-1       5-1   -1.355799190977935e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-1       5-2   -7.080820461285378e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-1       5-3    1.057221006713734e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-1       5-0   -1.549659949075030e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-1       5-1   -5.583811652924849e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-1       5-2   -1.996224487794656e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-1       5-3    4.774227217700482e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-1       5-0    8.864520385006748e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-1       5-1    1.701452334361741e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-1       5-2    1.092336916654610e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-1       5-3   -1.328567098916731e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-2       5-0   -5.850895297076977e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-2       5-1   -5.612905371621625e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-2       5-2   -2.109337056927748e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-2       5-3    6.227068117149350e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-2       5-0   -6.339931561521386e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-2       5-1   -7.080820461285364e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-2       5-2   -2.193037679558781e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-2       5-3    7.533896335956303e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-2       5-0   -2.662147958427581e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-2       5-1   -1.996224487794668e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-2       5-2   -9.508471542055009e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-2       5-3    2.216156257001325e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-2       5-0    7.958531960397842e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-2       5-1    1.092336916654588e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-2       5-2    2.754192107563037e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-2       5-3   -1.158260836205820e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-3       5-0   -3.163029162804833e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-3       5-1    1.071094448104043e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-3       5-2    6.227068117149306e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-0       5-3       5-3   -3.004508733275824e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-3       5-0   -6.070004790509594e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-3       5-1    1.057221006713734e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-3       5-2    7.533896335956983e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-1       5-3       5-3   -2.940388320173312e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-3       5-0   -1.119970267389317e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-3       5-1    4.774227217700482e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-3       5-2    2.216156257001309e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-2       5-3       5-3   -1.336716363690566e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-3       5-0    9.399142261962735e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-3       5-1   -1.328567098916731e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-3       5-2   -1.158260836205809e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-2       4-3       5-3       5-3    3.695315176045684e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-0       5-0    7.521537027004099e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-0       5-1    1.220827707549327e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-0       5-2    1.054749127222759e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-0       5-3    1.091629478003623e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-0       5-0    8.693532935491385e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-0       5-1    7.521484092497883e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-0       5-2    1.271104334664337e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-0       5-3    5.286945153289129e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-0       5-0    5.777752975611132e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-0       5-1    8.864520385006801e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-0       5-2    7.958531960397842e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-0       5-3    9.399142261968859e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-0       5-0   -2.782936849987281e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-0       5-1   -2.115448804617099e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-0       5-2   -4.038069821885032e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-0       5-3   -1.498638419724892e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-1       5-0    1.220827707549341e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-1       5-1    2.241395247038475e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-1       5-2    1.536318779421175e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-1       5-3   -1.799547610888852e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-1       5-0    7.521484092497820e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-1       5-1    2.664155641631492e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-1       5-2    9.714497984921449e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-1       5-3   -2.285369644867212e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-1       5-0    8.864520385006690e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-1       5-1    1.701452334361741e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-1       5-2    1.092336916654674e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-1       5-3   -1.328567098916731e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-1       5-0   -2.115448804617143e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-1       5-1   -8.484657043221032e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-1       5-2   -2.730374202614959e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-1       5-3    7.195484184223335e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-2       5-0    1.054749127222758e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-2       5-1    1.536318779421143e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-2       5-2    3.683955690278243e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-2       5-3   -1.657379041409383e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-2       5-0    1.271104334664337e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-2       5-1    9.714497984921558e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-2       5-2    4.545162992418934e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-2       5-3   -1.081392657106268e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-2       5-0    7.958531960397842e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-2       5-1    1.092336916654652e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-2       5-2    2.754192107563037e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-2       5-3   -1.158260836205809e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-2       5-0   -4.038069821885032e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-2       5-1   -2.730374202614970e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-2       5-2   -1.438122762846572e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-2       5-3    3.037548772980243e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-3       5-0    1.091629478003710e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-3       5-1   -1.799547610888851e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-3       5-2   -1.657379041409383e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-0       5-3       5-3    5.008884391539082e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-3       5-0    5.286945153289129e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-3       5-1   -2.285369644867212e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-3       5-2   -1.081392657106271e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-1       5-3       5-3    6.399605663311412e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-3       5-0    9.399142261958433e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-3       5-1   -1.328567098916731e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-3       5-2   -1.158260836205798e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-2       5-3       5-3    3.695315176045685e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-3       5-0   -1.498638419724413e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-3       5-1    7.195484184223337e+00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-3       5-2    3.037548772980216e-00\n";
    integralFileTwoBodyFADFingerPrint
        << "4-3       4-3       5-3       5-3   -2.013327838133568e+00\n";
    integralFileTwoBodyFADFingerPrint.close();
    //
    integralFileOneBodyWater.open("integral_file_test_OneBodyWater");
    integralFileOneBodyWater << "1-0    1-0    822.20658268\n";
    integralFileOneBodyWater << "1-1    1-1    2458.13980365\n";
    integralFileOneBodyWater << "1-2    1-2    4080.70800235\n";
    integralFileOneBodyWater << "1-3    1-3    5695.76754390\n";
    integralFileOneBodyWater << "1-4    1-4    7322.63669604\n";
    integralFileOneBodyWater << "1-5    1-5    9004.06198596\n";
    integralFileOneBodyWater << "2-0    2-0    1907.70566981\n";
    integralFileOneBodyWater << "2-1    2-1    5657.87641541\n";
    integralFileOneBodyWater << "2-2    2-2    9339.26135683\n";
    integralFileOneBodyWater << "2-3    2-3    13037.37641030\n";
    integralFileOneBodyWater << "2-4    2-4    16951.01097233\n";
    integralFileOneBodyWater << "2-5    2-5    21255.49635005\n";
    integralFileOneBodyWater << "3-0    3-0    1995.72280120\n";
    integralFileOneBodyWater << "3-1    3-1    6033.64526923\n";
    integralFileOneBodyWater << "3-2    3-2    10163.31053674\n";
    integralFileOneBodyWater << "3-3    3-3    14383.49653355\n";
    integralFileOneBodyWater << "3-4    3-4    18697.99018599\n";
    integralFileOneBodyWater << "3-5    3-5    23132.08260661\n";
    integralFileOneBodyWater.close();
    // Sets it now because integralsOneBodyFAD has to be populated before
    // passing it to the parameter container.
    parametersFADOneBodyBinary.set("L", 39);
    parametersFADOneBodyBinary.set("nmode_num_modes", 1);
    parametersFADOneBodyBinary.set("nmode_max_coupling", 1);
    parametersFADOneBodyBinary.set("nmode_num_basis", "39");
    parametersFADOneBodyBinary.set("symmetry", "nu1");
    parametersFADOneBodyBinary.set("LATTICE", "nmode lattice");
    parametersFADOneBodyBinary.set("MODEL", "nmode");
    parametersFADOneBodyBinary.set(
        "integrals_binary", maquis::serialize(integralsOneBodyFAD)
    );
    // Parameters for the fingerprint calculation
    parametersFADTwoBodyFingerPrint.set("L", 20);
    parametersFADTwoBodyFingerPrint.set("nmode_num_modes", 5);
    parametersFADTwoBodyFingerPrint.set("nmode_num_basis", "4,4,4,4,4");
    parametersFADTwoBodyFingerPrint.set("nmode_max_coupling", 2);
    parametersFADTwoBodyFingerPrint.set("symmetry", "nu1");
    parametersFADTwoBodyFingerPrint.set("LATTICE", "nmode lattice");
    parametersFADTwoBodyFingerPrint.set("MODEL", "nmode");
    parametersFADTwoBodyFingerPrint.set(
        "integral_file", "integral_file_test_TwoBodyFAD_Fingerprint"
    );
  }

  /** @brief Class destructor (removes tmp files) */
  ~NModeFixture() {
    std::remove("integral_file_test_OneBodyFAD");
    std::remove("integral_file_test_TwoBodyFAD");
    std::remove("integral_file_test_TwoBodyFAD_Fingerprint");
  }

  // Class members
  DmrgParameters parametersTwoMode, parametersFourMode, parametersFADOneBody, parametersFADOneBodyPaired,
      parametersFADTwoBody, parametersFADOneBodyBinary,
      parametersFADTwoBodyFingerPrint, parametersWater;
  std::ofstream integralFileOneBodyFAD, integralFileTwoBodyFADFingerPrint,
      integralFileTwoBodyFAD, integralFileOneBodyWater;
  MaquisIntegralType integralsOneBodyFAD;
};

#endif
