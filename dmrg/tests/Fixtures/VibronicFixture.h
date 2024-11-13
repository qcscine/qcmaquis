/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef TEST_VIBRONIC_FIXTURE_H
#define TEST_VIBRONIC_FIXTURE_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

/** @brief Fixture class for the test of the Vibronic DMRG */
struct VibronicFixture {
  /** @brief Constructor for the fixture class */
  VibronicFixture() {
    // Single-excited state excitonic model
    parametersExcitonicAggregate.set("L", 66);
    parametersExcitonicAggregate.set("symmetry", "u1");
    parametersExcitonicAggregate.set("LATTICE", "vibronic lattice");
    parametersExcitonicAggregate.set("MODEL", "excitonic");
    parametersExcitonicAggregate.set("Nmax", 8);
    parametersExcitonicAggregate.set("vibronic_num_elestates", 1);
    parametersExcitonicAggregate.set("vibronic_num_vibmodes", 10);
    parametersExcitonicAggregate.set("vibronic_num_molecules", 6);
    parametersExcitonicAggregate.set("vibronic_J_coupling", -500);
    parametersExcitonicAggregate.set("vibronic_J_interaction_type", "nn");
    parametersExcitonicAggregate.set(
        "integral_file", "integral_file_Excitonic"
    );
    parametersExcitonicAggregate.set("hamiltonian_units", "cm-1");
    //
    parametersExcitonicAggregateTwoSites.set("L", 22);
    parametersExcitonicAggregateTwoSites.set("symmetry", "u1");
    parametersExcitonicAggregateTwoSites.set("LATTICE", "vibronic lattice");
    parametersExcitonicAggregateTwoSites.set("MODEL", "excitonic");
    parametersExcitonicAggregateTwoSites.set("Nmax", 8);
    parametersExcitonicAggregateTwoSites.set("vibronic_num_elestates", 1);
    parametersExcitonicAggregateTwoSites.set("vibronic_num_vibmodes", 10);
    parametersExcitonicAggregateTwoSites.set("vibronic_num_molecules", 2);
    parametersExcitonicAggregateTwoSites.set("vibronic_J_coupling", -500);
    parametersExcitonicAggregateTwoSites.set(
        "vibronic_J_interaction_type", "nn"
    );
    parametersExcitonicAggregateTwoSites.set(
        "integral_file", "integral_file_Excitonic"
    );
    parametersExcitonicAggregateTwoSites.set("hamiltonian_units", "cm-1");
    // Vibronic models for S1/S2 states of pyrazine --> L=26
    parametersVibronic.set("L", 26);
    parametersVibronic.set("symmetry", "u1");
    parametersVibronic.set("LATTICE", "vibronic lattice");
    parametersVibronic.set("MODEL", "vibronic");
    parametersVibronic.set("Nmax", 6);
    parametersVibronic.set("vibronic_num_elestates", 2);
    parametersVibronic.set("vibronic_num_vibmodes", 24);
    parametersVibronic.set("hamiltonian_units", "cm-1");
    // "Fake" vibronic Hamiltonian with only a single state, with an Harmonic
    // PES.
    parametersFakeVibronic.set("L", 4);
    parametersFakeVibronic.set("symmetry", "u1");
    parametersFakeVibronic.set("LATTICE", "vibronic lattice");
    parametersFakeVibronic.set("MODEL", "vibronic");
    parametersFakeVibronic.set("Nmax", 6);
    parametersFakeVibronic.set("vibronic_num_elestates", 1);
    parametersFakeVibronic.set("vibronic_num_vibmodes", 3);
    parametersFakeVibronic.set("integral_file", "integral_file_VibronicFile");
    // Vibronic Hamiltonian for pyrazine, for the 4-mode Harmonic Hamiltonian.
    parametersVibronicPyrazineRedDim.set("L", 6);
    parametersVibronicPyrazineRedDim.set("symmetry", "u1");
    parametersVibronicPyrazineRedDim.set("LATTICE", "vibronic lattice");
    parametersVibronicPyrazineRedDim.set("MODEL", "vibronic");
    parametersVibronicPyrazineRedDim.set("Nmax", 6);
    parametersVibronicPyrazineRedDim.set("vibronic_num_elestates", 2);
    parametersVibronicPyrazineRedDim.set("vibronic_num_vibmodes", 4);
    parametersVibronicPyrazineRedDim.set(
        "integral_file", "integral_file_vibronic_Pyrazine_RedDim"
    );
    parametersVibronicPyrazineRedDim.set("hamiltonian_units", "cm-1");
    // Vibronic Hamiltonian for pyrazine, for the 4-mode full vibronic
    // Hamiltonian.
    parametersVibronicPyrazineRedDimFull.set("L", 6);
    parametersVibronicPyrazineRedDimFull.set("symmetry", "u1");
    parametersVibronicPyrazineRedDimFull.set("LATTICE", "vibronic lattice");
    parametersVibronicPyrazineRedDimFull.set("MODEL", "vibronic");
    parametersVibronicPyrazineRedDimFull.set("Nmax", 6);
    parametersVibronicPyrazineRedDimFull.set("vibronic_num_elestates", 2);
    parametersVibronicPyrazineRedDimFull.set("vibronic_num_vibmodes", 4);
    parametersVibronicPyrazineRedDimFull.set(
        "integral_file", "integral_file_vibronic_Pyrazine_RedDim_Full"
    );
    parametersVibronicPyrazineRedDimFull.set("hamiltonian_units", "cm-1");
    // Vibronic Hamiltonian for the thiphene dimer with a single harmonic mode
    parametersVibronicThiopheneDimer.set("L", 4);
    parametersVibronicThiopheneDimer.set("symmetry", "u1");
    parametersVibronicThiopheneDimer.set("LATTICE", "vibronic lattice");
    parametersVibronicThiopheneDimer.set("MODEL", "excitonic");
    parametersVibronicThiopheneDimer.set("Nmax", 6);
    parametersVibronicThiopheneDimer.set("vibronic_num_elestates", 1);
    parametersVibronicThiopheneDimer.set("vibronic_num_vibmodes", 1);
    parametersVibronicThiopheneDimer.set("vibronic_num_molecules", 2);
    parametersVibronicThiopheneDimer.set("vibronic_num_excitons", 1);
    parametersVibronicThiopheneDimer.set("vibronic_J_coupling", 0);
    parametersVibronicThiopheneDimer.set("vibronic_sorting", "intertwined");
    parametersVibronicThiopheneDimer.set(
        "integral_file", "integralFileThiopheneDimer"
    );
    parametersVibronicThiopheneDimer.set("nsweeps", 6);
    parametersVibronicThiopheneDimer.set("max_bond_dimension", 20);
    parametersVibronicThiopheneDimer.set("init_type", "basis_state_generic");
    parametersVibronicThiopheneDimer.set("init_basis_state", "1,0,0,0");
    parametersVibronicThiopheneDimer.set("simulation_type", "evolve");
    parametersVibronicThiopheneDimer.set("propagator_accuracy", 1.0E-10);
    parametersVibronicThiopheneDimer.set("propagator_maxiter", 10);
    parametersVibronicThiopheneDimer.set("time_step", 1);
    parametersVibronicThiopheneDimer.set("hamiltonian_units", "Hartree");
    parametersVibronicThiopheneDimer.set("time_units", "as");
    //
    parametersTestNmax.set("nsweeps", 10);
    parametersTestNmax.set("max_bond_dimension", 20);
    parametersTestNmax.set("optimization", "twosite");
    parametersTestNmax.set("integral_file", "integralFileTestNmax");
    parametersTestNmax.set("init_type", "basis_state_generic");
    parametersTestNmax.set("init_basis_state", "0,0,0,1,0,0,0,0");
    parametersTestNmax.set("Nmax", "1,2,3,4,5");
    parametersTestNmax.set("symmetry", "u1");
    parametersTestNmax.set("LATTICE", "vibronic lattice");
    parametersTestNmax.set("MODEL", "excitonicextended");
    parametersTestNmax.set("vibronic_J_coupling", 0);
    parametersTestNmax.set("vibronic_sorting", "intertwined");
    parametersTestNmax.set("vibronic_num_elestates", 1);
    parametersTestNmax.set("vibronic_num_vibmodes", 2);
    parametersTestNmax.set("vibronic_num_molecules", 3);
    parametersTestNmax.set("vibronic_num_excitons", 1);
    parametersTestNmax.set("vibronic_num_connectingmodes", 1);
    parametersTestNmax.set("L", 8);
    parametersTestNmax.set("simulation_type", "evolve");
    parametersTestNmax.set("propagator_accuracy", 1.0E-10);
    parametersTestNmax.set("propagator_maxiter", 10);
    parametersTestNmax.set("time_step", 1);
    parametersTestNmax.set("hamiltonian_units", "Hartree");
    parametersTestNmax.set("time_units", "as");
    //
    parametersExcitonicExtendedAggregate.set("L", 2);
    parametersExcitonicExtendedAggregate.set("symmetry", "u1");
    parametersExcitonicExtendedAggregate.set("LATTICE", "vibronic lattice");
    parametersExcitonicExtendedAggregate.set(
        "MODEL", "excitonic"
    );  // still have to change this...
    parametersExcitonicExtendedAggregate.set("Nmax", 8);
    parametersExcitonicExtendedAggregate.set("vibronic_num_elestates", 1);
    parametersExcitonicExtendedAggregate.set("vibronic_num_vibmodes", 1);
    parametersExcitonicExtendedAggregate.set("vibronic_num_molecules", 1);
    parametersExcitonicExtendedAggregate.set("vibronic_J_coupling", -2);
    parametersExcitonicExtendedAggregate.set(
        "integral_file", "integral_file_ExcitonicExtended"
    );
    parametersExcitonicExtendedAggregate.set("hamiltonian_units", "Hartree");
    //
    parametersSimpleCoherent.set("max_bond_dimension", 50);
    parametersSimpleCoherent.set(
        "integral_file", "integral_file_simpleCoherent"
    );
    parametersSimpleCoherent.set("L", 4);
    parametersSimpleCoherent.set("Nmax", 6);
    parametersSimpleCoherent.set("symmetry", "u1");
    parametersSimpleCoherent.set("LATTICE", "vibronic lattice");
    parametersSimpleCoherent.set("MODEL", "excitonicextended");
    parametersSimpleCoherent.set("vibronic_J_coupling", -0.0461);
    parametersSimpleCoherent.set("vibronic_sorting", "intertwined");
    parametersSimpleCoherent.set("vibronic_num_elestates", 1);
    parametersSimpleCoherent.set("vibronic_num_vibmodes", 1);
    parametersSimpleCoherent.set("vibronic_num_molecules", 2);
    parametersSimpleCoherent.set("vibronic_num_excitons", 1);
    parametersSimpleCoherent.set("vibronic_num_connectingmodes", 0);
    //
    parametersSimpleCoherentConnecting.set("nsweeps", 1);
    parametersSimpleCoherentConnecting.set("max_bond_dimension", 20);
    parametersSimpleCoherentConnecting.set("model_library", "coded");
    parametersSimpleCoherentConnecting.set("lattice_library", "coded");
    parametersSimpleCoherentConnecting.set("optimization", "singlesite");
    parametersSimpleCoherentConnecting.set(
        "integral_file", "integral_file_simpleCoherentConnecting"
    );
    parametersSimpleCoherentConnecting.set("Nmax", 6);
    parametersSimpleCoherentConnecting.set("symmetry", "u1");
    parametersSimpleCoherentConnecting.set("LATTICE", "vibronic lattice");
    parametersSimpleCoherentConnecting.set("MODEL", "excitonicextended");
    parametersSimpleCoherentConnecting.set("vibronic_j_coupling", 0.);
    parametersSimpleCoherentConnecting.set("vibronic_sorting", "intertwined");
    parametersSimpleCoherentConnecting.set("vibronic_num_elestates", 1);
    parametersSimpleCoherentConnecting.set("vibronic_num_vibmodes", 2);
    parametersSimpleCoherentConnecting.set("vibronic_num_molecules", 3);
    parametersSimpleCoherentConnecting.set("vibronic_num_excitons", 1);
    parametersSimpleCoherentConnecting.set("vibronic_num_connectingmodes", 1);
    parametersSimpleCoherentConnecting.set("L", 8);
    //
    parametersNmodeThiopheneOneBody.set("nsweeps", 1);
    parametersNmodeThiopheneOneBody.set("max_bond_dimension", 1);
    parametersNmodeThiopheneOneBody.set("model_library", "coded");
    parametersNmodeThiopheneOneBody.set("lattice_library", "coded");
    parametersNmodeThiopheneOneBody.set("optimization", "twosite");
    parametersNmodeThiopheneOneBody.set(
        "integral_file", "integral_file_NmodeThiopheneOneBody"
    );
    parametersNmodeThiopheneOneBody.set("Nmax", 6);
    parametersNmodeThiopheneOneBody.set("symmetry", "u1");
    parametersNmodeThiopheneOneBody.set("LATTICE", "vibronic lattice");
    parametersNmodeThiopheneOneBody.set("MODEL", "excitonicnmode");
    parametersNmodeThiopheneOneBody.set("vibronic_j_coupling", 0.);
    parametersNmodeThiopheneOneBody.set("vibronic_sorting", "intertwined");
    parametersNmodeThiopheneOneBody.set("vibronic_num_elestates", 1);
    parametersNmodeThiopheneOneBody.set("vibronic_num_vibmodes", 3);
    parametersNmodeThiopheneOneBody.set("vibronic_num_molecules", 10);
    parametersNmodeThiopheneOneBody.set("vibronic_num_excitons", 1);
    parametersNmodeThiopheneOneBody.set("vibronic_num_connectingmodes", 2);
    parametersNmodeThiopheneOneBody.set("L", 38);
    parametersNmodeThiopheneOneBody.set("init_type", "basis_state_generic");
    parametersNmodeThiopheneOneBody.set(
        "init_basis_state",
        "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
        "0,0,0"
    );
    parametersNmodeThiopheneOneBody.set("vibronic_max_coupling_nmode", 1);
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
    integralFilePyrazineRedDim << "-4114.23                0      0"
                               << std::endl;
    integralFilePyrazineRedDim << "  325.84799             1      1"
                               << std::endl;
    integralFilePyrazineRedDim << " -325.84799            -1     -1"
                               << std::endl;
    integralFilePyrazineRedDim << "  559.34550             2      2"
                               << std::endl;
    integralFilePyrazineRedDim << " -559.34550            -2     -2"
                               << std::endl;
    integralFilePyrazineRedDim << "  674.68277             3      3"
                               << std::endl;
    integralFilePyrazineRedDim << " -674.68277            -3     -3"
                               << std::endl;
    integralFilePyrazineRedDim << "  518.21124             4      4"
                               << std::endl;
    integralFilePyrazineRedDim << " -518.21124            -4     -4"
                               << std::endl;
    integralFilePyrazineRedDim << "EL_ST 1 1" << std::endl;
    integralFilePyrazineRedDim << " 4114.23                0      0"
                               << std::endl;
    integralFilePyrazineRedDim << "  325.84799             1      1"
                               << std::endl;
    integralFilePyrazineRedDim << " -325.84799            -1     -1"
                               << std::endl;
    integralFilePyrazineRedDim << "  559.34550             2      2"
                               << std::endl;
    integralFilePyrazineRedDim << " -559.34550            -2     -2"
                               << std::endl;
    integralFilePyrazineRedDim << "  674.68277             3      3"
                               << std::endl;
    integralFilePyrazineRedDim << " -674.68277            -3     -3"
                               << std::endl;
    integralFilePyrazineRedDim << "  518.21124             4      4"
                               << std::endl;
    integralFilePyrazineRedDim << " -518.21124            -4     -4"
                               << std::endl;
    //
    integralFilePyrazineRedDimFull.open(
        "integral_file_vibronic_Pyrazine_RedDim_Full"
    );
    integralFilePyrazineRedDimFull << "EL_ST 0 0" << std::endl;
    integralFilePyrazineRedDimFull << "-4114.23                0      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  325.84799             1      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -325.84799            -1     -1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  559.34550             2      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -559.34550            -2     -2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  674.68277             3      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -674.68277            -3     -3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  518.21124             4      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -518.21124            -4     -4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "EL_ST 1 1" << std::endl;
    integralFilePyrazineRedDimFull << " 4114.23                0      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  325.84799             1      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -325.84799            -1     -1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  559.34550             2      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -559.34550            -2     -2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  674.68277             3      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -674.68277            -3     -3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  518.21124             4      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -518.21124            -4     -4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "EL_ST 0 0" << std::endl;
    integralFilePyrazineRedDimFull << " -791.22989             1      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -405.69687             2      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " 1171.11703             3      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "EL_ST 1 1" << std::endl;
    integralFilePyrazineRedDimFull << " 1092.881251            1      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "-1379.877165            2      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  302.457910            3      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "EL_ST 0 0" << std::endl;
    integralFilePyrazineRedDimFull << "   0.16131088           1      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   8.71078783           1      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -16.45371035           1      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   8.71078783           2      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -65.33090875           2      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  38.23067993           2      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -16.45371035           3      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  38.23067993           3      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  -9.35603137           3      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -93.47965832           4      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "EL_ST 1 1 " << std::endl;
    integralFilePyrazineRedDimFull << " -73.96104114           1      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -24.03532198           1      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -15.24387871           1      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -24.03532198           2      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  39.35985614           2      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   9.27537593           2      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -15.24387871           3      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   9.27537593           3      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   1.77441974           3      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << " -93.47965832           4      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "EL_ST 0 1" << std::endl;
    integralFilePyrazineRedDimFull << " 1677.63321200          4      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  -80.65544294          4      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  -44.44114904          4      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   10.24324125          4      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   50.65161814          1      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  -44.44114904          2      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   10.24324125          3      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "EL_ST 1 0" << std::endl;
    integralFilePyrazineRedDimFull << " 1677.63321200          4      0"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  -80.65544294          4      1"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  -44.44114904          4      2"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   10.24324125          4      3"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   50.65161814          1      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "  -44.44114904          2      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull << "   10.24324125          3      4"
                                   << std::endl;
    integralFilePyrazineRedDimFull.close();
    //
    integralFileThiopheneDimer.open("integralFileThiopheneDimer");
    integralFileThiopheneDimer << "-0.002633     -1      -1" << std::endl;
    integralFileThiopheneDimer << "0.002633       1       1" << std::endl;
    integralFileThiopheneDimer << "0              1       0" << std::endl;
    integralFileThiopheneDimer.close();
    //
    integralFileExcitonicExtended.open("integral_file_ExcitonicExtended");
    integralFileExcitonicExtended << "0.1       0       0       -1      -1"
                                  << std::endl;
    integralFileExcitonicExtended << "0.2       0       0        1       1"
                                  << std::endl;
    integralFileExcitonicExtended << "0.3       1       0       -1      -1"
                                  << std::endl;
    integralFileExcitonicExtended << "0.4       1       0        0       0"
                                  << std::endl;
    integralFileExcitonicExtended << "0.5       1       0        1       0"
                                  << std::endl;
    integralFileExcitonicExtended << "0.6       1       0        1       1"
                                  << std::endl;
    integralFileExcitonic.close();
    //
    integralFileTestNmax.open("integralFileTestNmax");
    integralFileTestNmax << "-1. 0 0 -1 -1" << std::endl;
    integralFileTestNmax << "1. 0 0 1 1" << std::endl;
    integralFileTestNmax << "-1. 1 0 -1 -1" << std::endl;
    integralFileTestNmax << "1. 1 0 1 1" << std::endl;
    integralFileTestNmax << "-1.5 0 1 -2 -2" << std::endl;
    integralFileTestNmax << "1.5 0 1 2 2" << std::endl;
    integralFileTestNmax << "-2. 1 1 -2 -2" << std::endl;
    integralFileTestNmax << "2. 1 1 2 2" << std::endl;
    integralFileTestNmax.close();
    //
    IntegralFileSimpleCoherent.open("integral_file_simpleCoherent");
    IntegralFileSimpleCoherent << "-1.0 0 0 -1 -1" << std::endl;
    IntegralFileSimpleCoherent << "1.0 0 0 1 1" << std::endl;
    IntegralFileSimpleCoherent << "-2.0 1 0 -1 -1" << std::endl;
    IntegralFileSimpleCoherent << "2.0 1 0 1 1" << std::endl;
    IntegralFileSimpleCoherent.close();
    //
    IntegralFileSimpleCoherentConnecting.open(
        "integral_file_simpleCoherentConnecting"
    );
    IntegralFileSimpleCoherentConnecting << "-1. 0 0 -1 -1" << std::endl;
    IntegralFileSimpleCoherentConnecting << "-1. 1 0 -1 -1" << std::endl;
    IntegralFileSimpleCoherentConnecting << "1. 0 0 1 1" << std::endl;
    IntegralFileSimpleCoherentConnecting << "2. 1 0 1 1" << std::endl;
    IntegralFileSimpleCoherentConnecting << "-1. 0 1 -2 -2" << std::endl;
    IntegralFileSimpleCoherentConnecting << "-1. 1 1 -2 -2" << std::endl;
    IntegralFileSimpleCoherentConnecting << "1. 0 1 2 2" << std::endl;
    IntegralFileSimpleCoherentConnecting << "3. 1 1 2 2" << std::endl;
    IntegralFileSimpleCoherentConnecting.close();
    //
    integralFileNmodeThiopheneOneBody.open("integral_file_NmodeThiopheneOneBody"
    );
    integralFileNmodeThiopheneOneBody << "1-0 1-0 0 0 0.0024695226597438247"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-1 1-1 0 0 0.0074041535374397074"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-2 1-2 0 0 0.012332781742559487"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-3 1-3 0 0 0.01725540926093563"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-4 1-4 0 0 0.02217203807483142"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-5 1-5 0 0 0.02708267016294182"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-0 1-0 0 1 0.16924647028417006"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-1 1-1 0 1 0.1761473222364848"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-2 1-2 0 1 0.18279205365870574"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-3 1-3 0 1 0.18918066649230939"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-4 1-4 0 1 0.19531316255323344"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "1-5 1-5 0 1 0.20118954353188565"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-0 2-0 1 0 0.0007978743559850753"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-1 2-1 1 0 0.0023843970092003275"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-2 2-2 1 0 0.003958612739923099"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-3 2-3 1 0 0.00552052223246975"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-4 2-4 1 0 0.007070126162533611"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-5 2-5 1 0 0.008607425197185091"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-0 2-0 1 1 0.19183739943728595"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-1 2-1 1 1 0.19410099787740723"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-2 2-2 1 1 0.19634965522994724"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-3 2-3 1 1 0.19858337224515338"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-4 2-4 1 1 0.20080214966525123"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "2-5 2-5 1 1 0.20300598822444332"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-0 3-0 1 0 5.8623829606151595e-05"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-1 3-1 1 0 0.00018948539056790688"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-2 3-2 1 0 0.00034960590483998095"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-3 3-3 1 0 0.0004048709098722345"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-4 3-4 1 0 0.000525625275277624"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-5 3-5 1 0 0.0005279377622841062"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-0 3-0 1 1 0.1946748359451256"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-1 3-1 1 1 0.1954465394882818"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-2 3-2 1 1 0.1957090420358867"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-3 3-3 1 1 0.1962307810160739"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-4 3-4 1 1 0.19679482816770608"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody << "3-5 3-5 1 1 0.19702658069219867"
                                      << std::endl;
    integralFileNmodeThiopheneOneBody.close();
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
  DmrgParameters parametersVibronic, parametersFakeVibronic,
      parametersExcitonicAggregate, parametersExcitonicAggregateTwoSites,
      parametersVibronicPyrazineRedDim, parametersVibronicPyrazineRedDimFull,
      parametersVibronicThiopheneDimer, parametersExcitonicExtendedAggregate,
      parametersTestNmax, parametersSimpleCoherent,
      parametersSimpleCoherentConnecting, parametersNmodeThiopheneOneBody;
  std::ofstream integralFileFakeVibronic, integralFileExcitonic,
      integralFileExcitonicHarmonic, integralFilePyrazineRedDim,
      integralFilePyrazineRedDimFull, integralFileThiopheneDimer,
      integralFileExcitonicExtended, integralFileTestNmax,
      IntegralFileSimpleCoherent, IntegralFileSimpleCoherentConnecting,
      integralFileNmodeThiopheneOneBody;
};

#endif
