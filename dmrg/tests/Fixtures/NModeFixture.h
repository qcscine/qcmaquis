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
        // == PARAMETERS FOR THE ONE-MODE FAD CALCULATION ==
        // These data were generated based on the PES published by Bowman and based on 
        // a DVR basis set
        parametersFADOneBody.set("L", 39);
        parametersFADOneBody.set("nmode_num_modes", 1);
        parametersFADOneBody.set("nmode_max_coupling", 1);
        parametersFADOneBody.set("nmode_num_basis", "39");
        parametersFADOneBody.set("symmetry", "nu1");
        parametersFADOneBody.set("LATTICE", "nmode lattice");
        parametersFADOneBody.set("MODEL", "nmode");
        parametersFADOneBody.set("integral_file", "integral_file_test_FAD_OneBodyFAD");

        // == INPUT FILE CREATIONS ==
        integralFileOneBodyFAD.open("integral_file_test_FAD_OneBodyFAD");
        integralFileOneBodyFAD << "       1-0       1-0  -2.359242429009664e+03\n";
        integralFileOneBodyFAD << "       1-1       1-1  -2.358797386438359e+03\n";
        integralFileOneBodyFAD << "       1-2       1-2  -1.444541233850135e+03\n";
        integralFileOneBodyFAD << "       1-3       1-3  -1.437160009122911e+03\n";
        integralFileOneBodyFAD << "       1-4       1-4  -6.657064784181758e+03\n";
        integralFileOneBodyFAD << "       1-5       1-5  -6.128584661601011e+03\n";
        integralFileOneBodyFAD << "       1-6       1-6  -3.441604675890336e+03\n";
        integralFileOneBodyFAD << "       1-7       1-7   1.538026675462226e+03\n";
        integralFileOneBodyFAD << "       1-8       1-8   6.026905892315631e+03\n";
        integralFileOneBodyFAD << "       1-9       1-9   9.519079627438325e+03\n";
        integralFileOneBodyFAD << "      1-10      1-10   1.408367346423375e+03\n";
        integralFileOneBodyFAD << "      1-11      1-11   1.877962833530346e+03\n";
        integralFileOneBodyFAD << "      1-12      1-12   2.406532356803376e+03\n";
        integralFileOneBodyFAD << "      1-13      1-13   2.974933802551137e+03\n";
        integralFileOneBodyFAD << "      1-14      1-14   3.590607405365762e+03\n";
        integralFileOneBodyFAD << "      1-15      1-15   4.250019051366812e+03\n";
        integralFileOneBodyFAD << "      1-16      1-16   4.954168295956210e+03\n";
        integralFileOneBodyFAD << "      1-17      1-17   5.702284980641777e+03\n";
        integralFileOneBodyFAD << "      1-18      1-18   6.494407879568767e+03\n";
        integralFileOneBodyFAD << "      1-19      1-19   7.330063654159827e+03\n";
        integralFileOneBodyFAD << "      1-20      1-20   8.209556409406479e+03\n";
        integralFileOneBodyFAD << "      1-21      1-21   9.132461752263011e+03\n";
        integralFileOneBodyFAD << "      1-22      1-22   1.009901501528866e+03\n";
        integralFileOneBodyFAD << "      1-23      1-23   1.110869147300428e+03\n";
        integralFileOneBodyFAD << "      1-24      1-24   1.216206565203700e+03\n";
        integralFileOneBodyFAD << "      1-25      1-25   1.325854103259667e+03\n";
        integralFileOneBodyFAD << "      1-26      1-26   1.439864696031767e+03\n";
        integralFileOneBodyFAD << "      1-27      1-27   1.558157142548529e+03\n";
        integralFileOneBodyFAD << "      1-28      1-28   1.680827752856561e+03\n";
        integralFileOneBodyFAD << "      1-29      1-29   1.807804326901167e+03\n";
        integralFileOneBodyFAD << "      1-30      1-30   1.939221830929778e+03\n";
        integralFileOneBodyFAD << "      1-31      1-31   2.074605464868975e+03\n";
        integralFileOneBodyFAD << "      1-32      1-32   2.214735474735023e+03\n";
        integralFileOneBodyFAD << "      1-33      1-33   2.359908634757324e+03\n";
        integralFileOneBodyFAD << "      1-34      1-34   2.503011884645919e+03\n";
        integralFileOneBodyFAD << "      1-35      1-35   2.651668901080508e+03\n";
        integralFileOneBodyFAD << "      1-36      1-36   2.861283145084534e+03\n";
        integralFileOneBodyFAD << "      1-37      1-37   2.946541753544523e+03\n";
        integralFileOneBodyFAD << "      1-37      1-38   2.795969730154564e-13\n";
        integralFileOneBodyFAD << "      1-38       1-0  -4.721452208883015e-13\n";
        integralFileOneBodyFAD << "      1-38       1-1   9.167304342968547e-13\n";
        integralFileOneBodyFAD << "      1-38       1-2   3.627696363978469e-13\n";
        integralFileOneBodyFAD << "      1-38       1-3   2.368394571952736e-13\n";
        integralFileOneBodyFAD << "      1-38       1-4  -1.389061557427853e-13\n";
        integralFileOneBodyFAD << "      1-38       1-5   1.503577338928864e-13\n";
        integralFileOneBodyFAD << "      1-38       1-6  -3.152902036131743e-13\n";
        integralFileOneBodyFAD << "      1-38       1-7   2.076156246433921e-13\n";
        integralFileOneBodyFAD << "      1-38       1-8   2.296264501786515e-13\n";
        integralFileOneBodyFAD << "      1-38       1-9   8.566375343452281e-13\n";
        integralFileOneBodyFAD << "      1-38      1-10   2.391446450047095e-13\n";
        integralFileOneBodyFAD << "      1-38      1-11   2.260571271188797e-13\n";
        integralFileOneBodyFAD << "      1-38      1-12   6.781713813566390e-13\n";
        integralFileOneBodyFAD << "      1-38      1-13   2.974435883143153e-13\n";
        integralFileOneBodyFAD << "      1-38      1-14  -1.998820913472199e-13\n";
        integralFileOneBodyFAD << "      1-38      1-15   4.568733516507884e-13\n";
        integralFileOneBodyFAD << "      1-38      1-16   3.331368189120332e-13\n";
        integralFileOneBodyFAD << "      1-38      1-17  -7.138646119543568e-13\n";
        integralFileOneBodyFAD << "      1-38      1-18   2.569912603035684e-13\n";
        integralFileOneBodyFAD << "      1-38      1-19  -2.450935167709958e-13\n";
        integralFileOneBodyFAD << "      1-38      1-20   5.472962024983402e-13\n";
        integralFileOneBodyFAD << "      1-38      1-21   4.259392184660996e-13\n";
        integralFileOneBodyFAD << "      1-38      1-22   7.043464171282988e-13\n";
        integralFileOneBodyFAD << "      1-38      1-23   1.142183379126971e-13\n";
        integralFileOneBodyFAD << "      1-38      1-24  -7.923897192693361e-13\n";
        integralFileOneBodyFAD << "      1-38      1-25   5.663325921504564e-13\n";
        integralFileOneBodyFAD << "      1-38      1-26  -3.640709520967220e-13\n";
        integralFileOneBodyFAD << "      1-38      1-27  -1.832252504016183e-13\n";
        integralFileOneBodyFAD << "      1-38      1-28   1.903638965211618e-13\n";
        integralFileOneBodyFAD << "      1-38      1-29   2.141593835863070e-13\n";
        integralFileOneBodyFAD << "      1-38      1-30  -3.093413318468879e-13\n";
        integralFileOneBodyFAD << "      1-38      1-31  -6.377190533458921e-13\n";
        integralFileOneBodyFAD << "      1-38      1-32  -4.402165107051867e-13\n";
        integralFileOneBodyFAD << "      1-38      1-33  -3.854868904553527e-13\n";
        integralFileOneBodyFAD << "      1-38      1-34   8.685352778778008e-13\n";
        integralFileOneBodyFAD << "      1-38      1-35   1.182635707137718e-13\n";
        integralFileOneBodyFAD << "      1-38      1-36  -3.593118546836929e-13\n";
        integralFileOneBodyFAD << "      1-38      1-37   2.676992294828838e-13\n";
        integralFileOneBodyFAD << "      1-38      1-38   3.133897924626650e+03\n";
        integralFileOneBodyFAD.close();
    }

    /** @brief Class destructor (removes tmp files) */
    ~NModeFixture() {
        std::remove("integral_file_test_OneBody_FAD");
    }

    // Class members
    DmrgParameters parametersTwoMode, parametersFourMode, parametersFADOneBody;
    std::ofstream integralFileOneBodyFAD;
};

#endif