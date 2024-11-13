/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MODULE Transcorrelated

#include "alps/numeric/matrix/vector.hpp"
#include "dmrg/block_matrix/symmetry/2u1.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"
#include "utils/bindings.hpp"
#include <boost/filesystem/operations.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <string>
#include <vector>
#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "Fixtures/H2Fixture.h"
#include "Fixtures/TranscorrelatedFixture.h"

/**
 * @brief Checks consistency between conventional and transcorrelated
 * implementations.
 *
 * The check is done by giving as input to the transcorrelated Hamiltonian
 * constructor, a conventional Hamiltonian. Note that the FCIDUMP will be
 * different because tcDMRG does not assume eightfold symmetry, but only
 * twofold. So, also the MPO will be different and, as a consequence, it is not
 * trivial that the two calculations give the same energy.
 */
BOOST_FIXTURE_TEST_CASE(
    TestTCMolecular_H2_VersusConventional, TranscorrelatedFixture
) {
#if defined(HAVE_TwoU1) and defined(DMRG_TD)
  parametersH2Conventional_ConventionalFormat.set("nsweeps", 10);
  parametersH2Conventional_ConventionalFormat.set("max_bond_dimension", 100);
  maquis::DMRGInterface<double> interfaceConventional(
      parametersH2Conventional_ConventionalFormat
  );
  interfaceConventional.optimize();
  auto energy1 = interfaceConventional.energy();
  //
  parametersH2Conventional_TranscorrelatedFormat.set(
      "transcorrelated_nsweeps_TI", 5
  );
  parametersH2Conventional_TranscorrelatedFormat.set(
      "transcorrelated_nsweeps_TD", 5
  );
  parametersH2Conventional_TranscorrelatedFormat.set("max_bond_dimension", 100);
  parametersH2Conventional_TranscorrelatedFormat.set("time_step", 10.);
  parametersH2Conventional_TranscorrelatedFormat.set("propagator_maxiter", 10);
  parametersH2Conventional_TranscorrelatedFormat.set(
      "TD_backpropagation", "no"
  );
  parametersH2Conventional_TranscorrelatedFormat.set("time_units", "fs");
  maquis::DMRGInterface<double> interfaceTranscorrelated(
      parametersH2Conventional_TranscorrelatedFormat
  );
  interfaceTranscorrelated.runTranscorrelated();
  auto energy2 = interfaceConventional.energy();
  // The reference energy was generated with the UCISD module of PySCF
  BOOST_CHECK_CLOSE(energy1, -1.102429823850713, 1.0E-3);
  BOOST_CHECK_CLOSE(energy1, energy2, 1.0E-8);
#endif  // HAVE_TwoU1 and DMRG_TD
}

#if defined(HAVE_TwoU1) and defined(DMRG_TD)

/**
 * @brief Verify coherence between tcDMRG and FCI.
 * Here we take a conventional Hamiltonian, we encode it in the transcorrelated
 * format, and check that the energy is coherent with that of a "poor-man"
 * implementation of CISD.
 * Note that we do the check for both the conventional and the transcorrelated
 * format.
 */
BOOST_FIXTURE_TEST_CASE(
    TestTCMolecular_H2_VersusFullCI, TranscorrelatedFixture
) {
  int iCont = 0;
  parametersH2Conventional_TranscorrelatedFormat.set(
      "transcorrelated_nsweeps_TI", 5
  );
  parametersH2Conventional_TranscorrelatedFormat.set(
      "transcorrelated_nsweeps_TD", 5
  );
  parametersH2Conventional_TranscorrelatedFormat.set("max_bond_dimension", 100);
  parametersH2Conventional_TranscorrelatedFormat.set("time_step", 10.);
  parametersH2Conventional_TranscorrelatedFormat.set("propagator_maxiter", 10);
  parametersH2Conventional_TranscorrelatedFormat.set(
      "TD_backpropagation", "no"
  );
  parametersH2Conventional_TranscorrelatedFormat.set("time_units", "fs");
  // Does conventional TI and transcorrelated.
  for (auto& iParameter : std::vector<DmrgParameters>{
           parametersH2Conventional_ConventionalFormat,
           parametersH2Conventional_TranscorrelatedFormat
       }) {
    iParameter.set("nsweeps", 10);
    iParameter.set("max_bond_dimension", 100);
    maquis::DMRGInterface<double> interface(iParameter);
    if (iCont == 0) {
      interface.optimize();
    } else {
      interface.runTranscorrelated();
    }
    auto energyDMRG = interface.energy();
    // Hand-made Full-CI
    if (iCont == 1) {
      // Setting this keyword to "yes" enables constructing the MPO with the
      // transcorrelated code
      parametersH2Conventional_TranscorrelatedFormat.set(
          "transcorrelated_hamiltonian", "yes"
      );
      std::string integralFileName =
          parametersH2Conventional_TranscorrelatedFormat
              ["transcorrelated_integral_file"];
      parametersH2Conventional_TranscorrelatedFormat.set(
          "integral_file", integralFileName
      );
    }
    auto lattice = Lattice(iParameter);
    auto model = Model<matrix, TwoU1>(lattice, iParameter);
    auto mpo = make_mpo(lattice, model);
    std::vector<MPS<matrix, TwoU1>> vectorOfMPS;
    for (const auto& iString : fullCIDeterminantsH2) {
      iParameter.set("init_type", "hf");
      iParameter.set("hf_occ", iString);
      vectorOfMPS.emplace_back(
          lattice.size(), *(model.initializer(lattice, iParameter))
      );
    }
    matrix hamiltonianMatrix(vectorOfMPS.size(), vectorOfMPS.size(), 0.0);
    matrix eigenVectors(vectorOfMPS.size(), vectorOfMPS.size(), 0.0);
    alps::numeric::vector<double> eigenValues(vectorOfMPS.size(), 0.0);
    for (int iRow = 0; iRow < vectorOfMPS.size(); iRow++) {
      for (int iCol = 0; iCol < vectorOfMPS.size(); iCol++) {
        hamiltonianMatrix(iRow, iCol) =
            expval(vectorOfMPS[iRow], vectorOfMPS[iCol], mpo) /
            std::sqrt(norm(vectorOfMPS[iRow]) * norm(vectorOfMPS[iCol]));
      }
    }
    alps::numeric::syev(hamiltonianMatrix, eigenVectors, eigenValues);
    BOOST_CHECK_CLOSE(eigenValues[vectorOfMPS.size() - 1], energyDMRG, 1.0E-8);
    iCont += 1;
  }
}

/** @brief Same as above, but for the transcorrelated Hamiltonian */
BOOST_FIXTURE_TEST_CASE(
    TestTCMolecular_H2_VersusFullCI_Transcorrelated, TranscorrelatedFixture
) {
  parametersH2Transcorrelated.set("transcorrelated_nsweeps_TI", 10);
  parametersH2Transcorrelated.set("transcorrelated_nsweeps_TC", 30);
  parametersH2Transcorrelated.set("max_bond_dimension", 100);
  parametersH2Transcorrelated.set("time_step", 10.);
  parametersH2Transcorrelated.set("propagator_maxiter", 10);
  parametersH2Transcorrelated.set("TD_backpropagation", "no");
  parametersH2Transcorrelated.set("time_units", "fs");
  parametersH2Transcorrelated.set("init_type", "hf");
  parametersH2Transcorrelated.set("hf_occ", "4,1,1,1,1,1,1,1,1,1");
  parametersH2Transcorrelated.set("optimization", "twosite");
  maquis::DMRGInterface<double> interface(parametersH2Transcorrelated);
  interface.runTranscorrelated();
  auto energyDMRG = maquis::real(interface.energy());
  // Hand-made Full-CI
  parametersH2Transcorrelated.set("transcorrelated_hamiltonian", "yes");
  auto lattice = Lattice(parametersH2Transcorrelated);
  auto model = Model<matrix, TwoU1>(lattice, parametersH2Transcorrelated);
  auto mpo = make_mpo(lattice, model);
  std::vector<MPS<matrix, TwoU1>> vectorOfMPS;
  for (const auto& iString : fullCIDeterminantsH2) {
    parametersH2Transcorrelated.set("init_type", "hf");
    parametersH2Transcorrelated.set("hf_occ", iString);
    vectorOfMPS.emplace_back(
        lattice.size(),
        *(model.initializer(lattice, parametersH2Transcorrelated))
    );
  }
  cmatrix hamiltonianMatrix(vectorOfMPS.size(), vectorOfMPS.size(), 0.0);
  alps::numeric::vector<std::complex<double>> eigenValues(
      vectorOfMPS.size(), 0.0
  );
  for (int iRow = 0; iRow < vectorOfMPS.size(); iRow++) {
    for (int iCol = 0; iCol < vectorOfMPS.size(); iCol++) {
      hamiltonianMatrix(iRow, iCol) =
          expval(vectorOfMPS[iRow], vectorOfMPS[iCol], mpo) /
          std::sqrt(norm(vectorOfMPS[iRow]) * norm(vectorOfMPS[iCol]));
    }
  }
  alps::numeric::geev(hamiltonianMatrix, eigenValues);
  // Checks that the eigenvalues are real (this comes from the fact that the
  // matrix is obtained as similarity transformation of a real-valued matrix)
  double minimumEnergy = std::real(eigenValues[0]);
  for (int iElement = 0; iElement < vectorOfMPS.size(); iElement++) {
    BOOST_CHECK_SMALL(std::imag(eigenValues[iElement]), 1.0E-8);
    if (iElement != 0) {
      auto realEnergy = std::real(eigenValues[iElement]);
      if (realEnergy < minimumEnergy) {
        minimumEnergy = realEnergy;
      }
    }
  }
  BOOST_CHECK_CLOSE(minimumEnergy, energyDMRG, 1.0E-8);
}

/**
 * @brief Transcorrelated DMRG calculation on Be.
 * The reference energy has been generated in this case with the Owl CC
 * code by Max Moerchen.
 */
BOOST_FIXTURE_TEST_CASE(TestTCMolecular_Be_VersusCC, TranscorrelatedFixture) {
  parametersBeTranscorrelatedTwoBody.set("propagator_accuracy", 1.0E-10);
  parametersBeTranscorrelatedTwoBody.set("propagator_maxiter", 10);
  parametersBeTranscorrelatedTwoBody.set("hamiltonian_units", "Hartree");
  parametersBeTranscorrelatedTwoBody.set("time_units", "fs");
  parametersBeTranscorrelatedTwoBody.set("TD_backpropagation", "no");
  parametersBeTranscorrelatedTwoBody.set("symmetry", "2u1");
  parametersBeTranscorrelatedTwoBody.set("chkpfile", "Be.tcDMRG.checkpoint.h5");
  parametersBeTranscorrelatedTwoBody.set("transcorrelated_nsweeps_TI", 0);
  parametersBeTranscorrelatedTwoBody.set(
      "integral_file", "IntegralFile_Be_Conventional"
  );
  parametersBeTranscorrelatedTwoBody.set(
      "transcorrelated_integral_file", "IntegralFile_Be_Transcorrelated_TwoBody"
  );

  parametersBeTranscorrelatedTwoBody.set("transcorrelated_nsweeps_TC", 10);
  parametersBeTranscorrelatedTwoBody.set("time_step", 0.1);
  maquis::DMRGInterface<double> interface(parametersBeTranscorrelatedTwoBody);
  interface.runTranscorrelated();
  parametersBeTranscorrelatedTwoBody.set("transcorrelated_nsweeps_TC", 20);
  parametersBeTranscorrelatedTwoBody.set("time_step", 0.01);
  // Reloads from checkpoint
  maquis::DMRGInterface<double> interface2(parametersBeTranscorrelatedTwoBody);
  interface2.runTranscorrelated();
  parametersBeTranscorrelatedTwoBody.set("transcorrelated_nsweeps_TC", 3);
  parametersBeTranscorrelatedTwoBody.set("time_step", 0.001);
  maquis::DMRGInterface<double> interface3(parametersBeTranscorrelatedTwoBody);
  interface3.runTranscorrelated();
  boost::filesystem::remove_all("Be.tcDMRG.checkpoint.h5");
  auto energyDMRG = maquis::real(interface3.energy());
  auto refEnergy = -14.6505807967243;
  BOOST_CHECK_SMALL(std::abs(energyDMRG - refEnergy), 1.0E-10);
}

/** @brief Checks quantum format for Hermitian Hamiltonians */
BOOST_FIXTURE_TEST_CASE(TestTCMolecular_H2_QuantumFormat, H2Fixture) {
  // Types definition
  using ModelType = Model<matrix, TwoU1>;
  using MPSType = MPS<matrix, TwoU1>;
  // Generates the conventional MPO
  auto lattice = Lattice(parametersH2QuantumFormatTranscorrelated);
  auto conventionalModel =
      ModelType(lattice, parametersH2QuantumFormatTranscorrelated);
  auto conventionalMpo = make_mpo(lattice, conventionalModel);
  parametersH2QuantumFormatTranscorrelated.set("init_type", "default");
  // Generates the MPS
  auto mps = MPSType(
      lattice.size(), *(conventionalModel.initializer(
                          lattice, parametersH2QuantumFormatTranscorrelated
                      ))
  );
  // Generates the transcorrelatedMPO
  auto transcorrelatedParametersContainer =
      parametersH2QuantumFormatTranscorrelated;
  transcorrelatedParametersContainer.set("transcorrelated_hamiltonian", "yes");
  transcorrelatedParametersContainer.set("imaginary_time", "yes");
  auto transcorrelatedModel =
      ModelType(lattice, transcorrelatedParametersContainer);
  auto transcorrelatedMpo = make_mpo(lattice, transcorrelatedModel);
  // Compares the energy calculated based on the two Hamiltonians
  auto energyConventional = expval(mps, conventionalMpo);
  auto energyQuantum = expval(mps, transcorrelatedMpo);
  BOOST_CHECK_CLOSE(energyConventional, energyQuantum, 1.0E-13);
  // Now also tries conventional Hamiltonian in quantum format
  auto conventionalModelQuantumFormat =
      ModelType(lattice, parametersH2QuantumFormat);
  auto conventionalMPOQuantumFormat =
      make_mpo(lattice, conventionalModelQuantumFormat);
  auto energyConventionalQuantumFormat =
      expval(mps, conventionalMPOQuantumFormat);
  BOOST_CHECK_CLOSE(energyConventionalQuantumFormat, energyQuantum, 1.0E-13);
}

/** @brief Checks quantum format for non-Hermitian Hamiltonians */
BOOST_FIXTURE_TEST_CASE(
    TestTCMolecular_H2_QuantumFormat_NonHermitian, TranscorrelatedFixture
) {
  // Types definition
  using ModelType = Model<matrix, TwoU1>;
  using MPSType = MPS<matrix, TwoU1>;
  // Generates the MPS
  parametersH2Transcorrelated.set("init_type", "hf");
  parametersH2Transcorrelated.set("hf_occ", "4,1,1,1,1,1,1,1,1,1");
  auto lattice = Lattice(parametersH2Transcorrelated);
  auto modelForMPS = ModelType(lattice, parametersH2Transcorrelated);
  auto randomMps = MPSType(
      lattice.size(),
      *(modelForMPS.initializer(lattice, parametersH2Transcorrelated))
  );
  // Generates the transcorrelatedMPO with the conventional format
  auto transcorrelatedParametersContainer = parametersH2Transcorrelated;
  transcorrelatedParametersContainer.set("transcorrelated_hamiltonian", "yes");
  transcorrelatedParametersContainer.set("imaginary_time", "yes");
  auto transcorrelatedModel =
      ModelType(lattice, transcorrelatedParametersContainer);
  auto transcorrelatedMpo = make_mpo(lattice, transcorrelatedModel);
  // Generates the transcorrelatedMPO in the quantum format
  transcorrelatedParametersContainer = parametersH2TranscorrelatedQuantumFormat;
  transcorrelatedParametersContainer.set("transcorrelated_hamiltonian", "yes");
  transcorrelatedParametersContainer.set("imaginary_time", "yes");
  auto transcorrelatedModelQF =
      ModelType(lattice, transcorrelatedParametersContainer);
  auto transcorrelatedMpoQF = make_mpo(lattice, transcorrelatedModelQF);
  // Checks coherence in the energy
  auto energyConventionalFormat = expval(randomMps, transcorrelatedMpo);
  auto energyQuantumFormat = expval(randomMps, transcorrelatedMpoQF);
  BOOST_CHECK_CLOSE(energyConventionalFormat, energyQuantumFormat, 1.0E-13);
}

#endif  // HAVE_TwoU1 and DMRG_TD
