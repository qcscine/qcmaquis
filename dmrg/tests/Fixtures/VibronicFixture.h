/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef TEST_VIBRONIC_FIXTURE_H
#define TEST_VIBRONIC_FIXTURE_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

/** @brief Fixture class for the test of the Vibronic DMRG */
struct VibronicFixture
{
    /** @brief Constructor for the fixture class */
    VibronicFixture() {
        // Single-excited state excitonic model
        parametersExcitonicAggregate.set("L", 66);
        parametersExcitonicAggregate.set("symmetry", "u1");
        parametersExcitonicAggregate.set("LATTICE", "vibronic lattice");
        parametersExcitonicAggregate.set("MODEL", "excitonic");
        parametersExcitonicAggregate.set("Nmax", 8);
        parametersExcitonicAggregate.set("vibronic_nstates", 1);
        parametersExcitonicAggregate.set("vibronic_nmodes", 10);
        parametersExcitonicAggregate.set("n_excitons", 6);
        parametersExcitonicAggregate.set("J_coupling", -500);
        parametersExcitonicAggregate.set("J_interaction", "nn");
        parametersExcitonicAggregate.set("integral_file", "integral_file_Excitonic");
        //
        parametersExcitonicAggregateTwoSites.set("L", 22);
        parametersExcitonicAggregateTwoSites.set("symmetry", "u1");
        parametersExcitonicAggregateTwoSites.set("LATTICE", "vibronic lattice");
        parametersExcitonicAggregateTwoSites.set("MODEL", "excitonic");
        parametersExcitonicAggregateTwoSites.set("Nmax", 8);
        parametersExcitonicAggregateTwoSites.set("vibronic_nstates", 1);
        parametersExcitonicAggregateTwoSites.set("vibronic_nmodes", 10);
        parametersExcitonicAggregateTwoSites.set("n_excitons", 2);
        parametersExcitonicAggregateTwoSites.set("J_coupling", -500);
        parametersExcitonicAggregateTwoSites.set("J_interaction", "nn");
        parametersExcitonicAggregateTwoSites.set("integral_file", "integral_file_Excitonic");
        // Vibronic models for S1/S2 states of pyrazine --> L=26
        parametersVibronic.set("L", 26);
        parametersVibronic.set("symmetry", "u1");
        parametersVibronic.set("LATTICE", "vibronic lattice");
        parametersVibronic.set("MODEL", "vibronic");
        parametersVibronic.set("Nmax", 6);
        parametersVibronic.set("vibronic_nstates", 2);
        parametersVibronic.set("vibronic_nmodes", 24);
        // "Fake" vibronic Hamiltonian with only a single state, with an Harmonic PES.
        parametersFakeVibronic.set("L", 4);
        parametersFakeVibronic.set("symmetry", "u1");
        parametersFakeVibronic.set("LATTICE", "vibronic lattice");
        parametersFakeVibronic.set("MODEL", "vibronic");
        parametersFakeVibronic.set("Nmax", 6);
        parametersFakeVibronic.set("vibronic_nstates", 1);
        parametersFakeVibronic.set("vibronic_nmodes", 3);
        parametersFakeVibronic.set("integral_file", "integral_file_VibronicFile");
        // Vibronic Hamiltonian for pyrazine, for the 4-mode Harmonic Hamiltonian.
        parametersVibronicPyrazineRedDim.set("L", 6);
        parametersVibronicPyrazineRedDim.set("symmetry", "u1");
        parametersVibronicPyrazineRedDim.set("LATTICE", "vibronic lattice");
        parametersVibronicPyrazineRedDim.set("MODEL", "vibronic");
        parametersVibronicPyrazineRedDim.set("Nmax", 6);
        parametersVibronicPyrazineRedDim.set("vibronic_nstates", 2);
        parametersVibronicPyrazineRedDim.set("vibronic_nmodes", 4);
        parametersVibronicPyrazineRedDim.set("integral_file", "integral_file_vibronic_Pyrazine_RedDim");
        // Vibronic Hamiltonian for pyrazine, for the 4-mode full vibronic Hamiltonian.
        parametersVibronicPyrazineRedDimFull.set("L", 6);
        parametersVibronicPyrazineRedDimFull.set("symmetry", "u1");
        parametersVibronicPyrazineRedDimFull.set("LATTICE", "vibronic lattice");
        parametersVibronicPyrazineRedDimFull.set("MODEL", "vibronic");
        parametersVibronicPyrazineRedDimFull.set("Nmax", 6);
        parametersVibronicPyrazineRedDimFull.set("vibronic_nstates", 2);
        parametersVibronicPyrazineRedDimFull.set("vibronic_nmodes", 4);
        parametersVibronicPyrazineRedDimFull.set("integral_file", "integral_file_vibronic_Pyrazine_RedDim_Full");
        //
        integralFileFakeVibronic.open("integral_file_VibronicFile");
        integralFileFakeVibronic << "EL_ST 0 0 " << std::endl;
        integralFileFakeVibronic << " 1.0000   1   1  " << std::endl;
        integralFileFakeVibronic << "-1.0000  -1  -1  " << std::endl;
        integralFileFakeVibronic << " 2.0000   2   2  " << std::endl;
        integralFileFakeVibronic << "-2.0000  -2  -2  " << std::endl;
        integralFileFakeVibronic << " 3.0000   3   3  " << std::endl;
        integralFileFakeVibronic << "-3.0000  -3  -3  " << std::endl;
        integralFileFakeVibronic.close();
        //
        integralFileExcitonic.open("integral_file_Excitonic");
        integralFileExcitonic << " 103.00      1     1" << std::endl;
        integralFileExcitonic << "-103.00     -1    -1" << std::endl;
        integralFileExcitonic << " 105.50      2     2" << std::endl;
        integralFileExcitonic << "-105.50     -2    -2" << std::endl;
        integralFileExcitonic << " 270.00      3     3" << std::endl;
        integralFileExcitonic << "-270.00     -3    -3" << std::endl;
        integralFileExcitonic << " 276.00      4     4" << std::endl;
        integralFileExcitonic << "-276.00     -4    -4" << std::endl;
        integralFileExcitonic << " 375.50      5     5" << std::endl;
        integralFileExcitonic << "-375.50     -5    -5" << std::endl;
        integralFileExcitonic << " 662.50      6     6" << std::endl;
        integralFileExcitonic << "-662.50     -6    -6" << std::endl;
        integralFileExcitonic << " 685.50      7     7" << std::endl;
        integralFileExcitonic << "-685.50     -7    -7" << std::endl;
        integralFileExcitonic << " 734.50      8     8" << std::endl;
        integralFileExcitonic << "-734.50     -8    -8" << std::endl;
        integralFileExcitonic << " 785.50      9     9" << std::endl;
        integralFileExcitonic << "-785.50     -9    -9" << std::endl;
        integralFileExcitonic << " 814.50     10    10" << std::endl;
        integralFileExcitonic << "-814.50    -10   -10" << std::endl;
        integralFileExcitonic << " 129.30      1     0" << std::endl;
        integralFileExcitonic << " 138.36      2     0" << std::endl;
        integralFileExcitonic << " 105.26      3     0" << std::endl;
        integralFileExcitonic << " 150.16      4     0" << std::endl;
        integralFileExcitonic << " 192.93      5     0" << std::endl;
        integralFileExcitonic << " 187.38      6     0" << std::endl;
        integralFileExcitonic << " 884.27      7     0" << std::endl;
        integralFileExcitonic << " 425.76      8     0" << std::endl;
        integralFileExcitonic << " 639.67      9     0" << std::endl;
        integralFileExcitonic << " 454.68     10     0" << std::endl;
        integralFileExcitonic.close();
        //
        integralFileExcitonicHarmonic.open("integral_file_Excitonic_Harmonic");
        integralFileExcitonicHarmonic << " 103.00      1     1" << std::endl;
        integralFileExcitonicHarmonic << "-103.00     -1    -1" << std::endl;
        integralFileExcitonicHarmonic << " 105.50      2     2" << std::endl;
        integralFileExcitonicHarmonic << "-105.50     -2    -2" << std::endl;
        integralFileExcitonicHarmonic << " 270.00      3     3" << std::endl;
        integralFileExcitonicHarmonic << "-270.00     -3    -3" << std::endl;
        integralFileExcitonicHarmonic << " 276.00      4     4" << std::endl;
        integralFileExcitonicHarmonic << "-276.00     -4    -4" << std::endl;
        integralFileExcitonicHarmonic << " 375.50      5     5" << std::endl;
        integralFileExcitonicHarmonic << "-375.50     -5    -5" << std::endl;
        integralFileExcitonicHarmonic << " 662.50      6     6" << std::endl;
        integralFileExcitonicHarmonic << "-662.50     -6    -6" << std::endl;
        integralFileExcitonicHarmonic << " 685.50      7     7" << std::endl;
        integralFileExcitonicHarmonic << "-685.50     -7    -7" << std::endl;
        integralFileExcitonicHarmonic << " 734.50      8     8" << std::endl;
        integralFileExcitonicHarmonic << "-734.50     -8    -8" << std::endl;
        integralFileExcitonicHarmonic << " 785.50      9     9" << std::endl;
        integralFileExcitonicHarmonic << "-785.50     -9    -9" << std::endl;
        integralFileExcitonicHarmonic << " 814.50     10    10" << std::endl;
        integralFileExcitonicHarmonic << "-814.50    -10   -10" << std::endl;
        integralFileExcitonicHarmonic.close();
        //
        integralFilePyrazineRedDim.open("integral_file_vibronic_Pyrazine_RedDim");
        integralFilePyrazineRedDim << "EL_ST 0 0" << std::endl;
        integralFilePyrazineRedDim << "-4114.23                0      0" << std::endl;
        integralFilePyrazineRedDim << "  325.84799             1      1" << std::endl;
        integralFilePyrazineRedDim << " -325.84799            -1     -1" << std::endl;
        integralFilePyrazineRedDim << "  559.34550             2      2" << std::endl;
        integralFilePyrazineRedDim << " -559.34550            -2     -2" << std::endl;
        integralFilePyrazineRedDim << "  674.68277             3      3" << std::endl;
        integralFilePyrazineRedDim << " -674.68277            -3     -3" << std::endl;
        integralFilePyrazineRedDim << "  518.21124             4      4" << std::endl;
        integralFilePyrazineRedDim << " -518.21124            -4     -4" << std::endl;
        integralFilePyrazineRedDim << "EL_ST 1 1" << std::endl;
        integralFilePyrazineRedDim << " 4114.23                0      0" << std::endl;
        integralFilePyrazineRedDim << "  325.84799             1      1" << std::endl;
        integralFilePyrazineRedDim << " -325.84799            -1     -1" << std::endl;
        integralFilePyrazineRedDim << "  559.34550             2      2" << std::endl;
        integralFilePyrazineRedDim << " -559.34550            -2     -2" << std::endl;
        integralFilePyrazineRedDim << "  674.68277             3      3" << std::endl;
        integralFilePyrazineRedDim << " -674.68277            -3     -3" << std::endl;
        integralFilePyrazineRedDim << "  518.21124             4      4" << std::endl;
        integralFilePyrazineRedDim << " -518.21124            -4     -4" << std::endl;
        //
        integralFilePyrazineRedDimFull.open("integral_file_vibronic_Pyrazine_RedDim_Full");
        integralFilePyrazineRedDimFull << "EL_ST 0 0" << std::endl;
        integralFilePyrazineRedDimFull << "-4114.23                0      0" << std::endl;
        integralFilePyrazineRedDimFull << "  325.84799             1      1" << std::endl;
        integralFilePyrazineRedDimFull << " -325.84799            -1     -1" << std::endl;
        integralFilePyrazineRedDimFull << "  559.34550             2      2" << std::endl;
        integralFilePyrazineRedDimFull << " -559.34550            -2     -2" << std::endl;
        integralFilePyrazineRedDimFull << "  674.68277             3      3" << std::endl;
        integralFilePyrazineRedDimFull << " -674.68277            -3     -3" << std::endl;
        integralFilePyrazineRedDimFull << "  518.21124             4      4" << std::endl;
        integralFilePyrazineRedDimFull << " -518.21124            -4     -4" << std::endl;
        integralFilePyrazineRedDimFull << "EL_ST 1 1" << std::endl;
        integralFilePyrazineRedDimFull << " 4114.23                0      0" << std::endl;
        integralFilePyrazineRedDimFull << "  325.84799             1      1" << std::endl;
        integralFilePyrazineRedDimFull << " -325.84799            -1     -1" << std::endl;
        integralFilePyrazineRedDimFull << "  559.34550             2      2" << std::endl;
        integralFilePyrazineRedDimFull << " -559.34550            -2     -2" << std::endl;
        integralFilePyrazineRedDimFull << "  674.68277             3      3" << std::endl;
        integralFilePyrazineRedDimFull << " -674.68277            -3     -3" << std::endl;
        integralFilePyrazineRedDimFull << "  518.21124             4      4" << std::endl;
        integralFilePyrazineRedDimFull << " -518.21124            -4     -4" << std::endl;
        integralFilePyrazineRedDimFull << "EL_ST 0 0" << std::endl;
        integralFilePyrazineRedDimFull << " -791.22989             1      0" << std::endl;
        integralFilePyrazineRedDimFull << " -405.69687             2      0" << std::endl;
        integralFilePyrazineRedDimFull << " 1171.11703             3      0" << std::endl;
        integralFilePyrazineRedDimFull << "EL_ST 1 1" << std::endl;
        integralFilePyrazineRedDimFull << " 1092.881251            1      0" << std::endl;
        integralFilePyrazineRedDimFull << "-1379.877165            2      0" << std::endl;
        integralFilePyrazineRedDimFull << "  302.457910            3      0" << std::endl;
        integralFilePyrazineRedDimFull << "EL_ST 0 0" << std::endl;
        integralFilePyrazineRedDimFull << "   0.16131088           1      1" << std::endl;
        integralFilePyrazineRedDimFull << "   8.71078783           1      2" << std::endl;
        integralFilePyrazineRedDimFull << " -16.45371035           1      3" << std::endl;
        integralFilePyrazineRedDimFull << "   8.71078783           2      1" << std::endl;
        integralFilePyrazineRedDimFull << " -65.33090875           2      2" << std::endl;
        integralFilePyrazineRedDimFull << "  38.23067993           2      3" << std::endl;
        integralFilePyrazineRedDimFull << " -16.45371035           3      1" << std::endl;
        integralFilePyrazineRedDimFull << "  38.23067993           3      2" << std::endl;
        integralFilePyrazineRedDimFull << "  -9.35603137           3      3" << std::endl;
        integralFilePyrazineRedDimFull << " -93.47965832           4      4" << std::endl;
        integralFilePyrazineRedDimFull << "EL_ST 1 1 " << std::endl;
        integralFilePyrazineRedDimFull << " -73.96104114           1      1" << std::endl;
        integralFilePyrazineRedDimFull << " -24.03532198           1      2" << std::endl;
        integralFilePyrazineRedDimFull << " -15.24387871           1      3" << std::endl;
        integralFilePyrazineRedDimFull << " -24.03532198           2      1" << std::endl;
        integralFilePyrazineRedDimFull << "  39.35985614           2      2" << std::endl;
        integralFilePyrazineRedDimFull << "   9.27537593           2      3" << std::endl;
        integralFilePyrazineRedDimFull << " -15.24387871           3      1" << std::endl;
        integralFilePyrazineRedDimFull << "   9.27537593           3      2" << std::endl;
        integralFilePyrazineRedDimFull << "   1.77441974           3      3" << std::endl;
        integralFilePyrazineRedDimFull << " -93.47965832           4      4" << std::endl;
        integralFilePyrazineRedDimFull << "EL_ST 0 1" << std::endl;
        integralFilePyrazineRedDimFull << " 1677.63321200          4      0" << std::endl;
        integralFilePyrazineRedDimFull << "  -80.65544294          4      1" << std::endl;
        integralFilePyrazineRedDimFull << "  -44.44114904          4      2" << std::endl;
        integralFilePyrazineRedDimFull << "   10.24324125          4      3" << std::endl;
        integralFilePyrazineRedDimFull << "   50.65161814          1      4" << std::endl;
        integralFilePyrazineRedDimFull << "  -44.44114904          2      4" << std::endl;
        integralFilePyrazineRedDimFull << "   10.24324125          3      4" << std::endl;
        integralFilePyrazineRedDimFull << "EL_ST 1 0" << std::endl;
        integralFilePyrazineRedDimFull << " 1677.63321200          4      0" << std::endl;
        integralFilePyrazineRedDimFull << "  -80.65544294          4      1" << std::endl;
        integralFilePyrazineRedDimFull << "  -44.44114904          4      2" << std::endl;
        integralFilePyrazineRedDimFull << "   10.24324125          4      3" << std::endl;
        integralFilePyrazineRedDimFull << "   50.65161814          1      4" << std::endl;
        integralFilePyrazineRedDimFull << "  -44.44114904          2      4" << std::endl;
        integralFilePyrazineRedDimFull << "   10.24324125          3      4" << std::endl;
        integralFilePyrazineRedDimFull.close();
    }

    /** @brief Class destructor */
    ~VibronicFixture() {
        std::remove("integral_file_VibronicFile");
        std::remove("integral_file_Excitonic");
        std::remove("integral_file_Excitonic_Harmonic");
        std::remove("integral_file_vibronic_Pyrazine_RedDim");
        std::remove("integral_file_vibronic_Pyrazine_RedDim_Full");
    }

    // Class members
    DmrgParameters parametersVibronic, parametersFakeVibronic, parametersExcitonicAggregate,
        parametersExcitonicAggregateTwoSites, parametersVibronicPyrazineRedDim, parametersVibronicPyrazineRedDimFull;
    std::ofstream integralFileFakeVibronic, integralFileExcitonic, integralFileExcitonicHarmonic,
        integralFilePyrazineRedDim, integralFilePyrazineRedDimFull;
};

#endif
