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

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include <boost/test/included/unit_test.hpp>

#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"

#include "Fixtures/TranscorrelatedFixture.h"

/**
 * @brief Checks consistency between conventional and transcorrelated implementations.
 * 
 * The check is done by giving as input to the transcorrelated Hamiltonian constructor,
 * a conventional Hamiltonian. Note that the FCIDUMP will be different because tcDMRG
 * does not assume eightfold symmetry, but only twofold. So, also the MPO will be different
 * and, as a consequence, it is not trivial that the two calculations give the same energy.
 */
BOOST_FIXTURE_TEST_CASE(TestTCMolecular_H2_VersusConventional, TranscorrelatedFixture)
{
#ifdef HAVE_TwoU1
    parametersH2Conventional_ConventionalFormat.set("nsweeps", 10);
    parametersH2Conventional_ConventionalFormat.set("max_bond_dimension", 100);
    maquis::DMRGInterface<double> interfaceConventional(parametersH2Conventional_ConventionalFormat);
    interfaceConventional.optimize();
    auto energy1 = interfaceConventional.energy();
    //
    parametersH2Conventional_TranscorrelatedFormat.set("nsweeps", 10);
    parametersH2Conventional_TranscorrelatedFormat.set("max_bond_dimension", 100);
    maquis::DMRGInterface<double> interfaceTranscorrelated(parametersH2Conventional_TranscorrelatedFormat);
    interfaceTranscorrelated.optimize();
    auto energy2 = interfaceConventional.energy();
    // The reference energy was generated with the UCISD module of PySCF
    BOOST_CHECK_CLOSE(energy1, -1.102429823850713, 1.0E-3);
    BOOST_CHECK_CLOSE(energy1, energy2, 1.0E-8);
#endif
}

#ifdef HAVE_TwoU1

/**
 * @brief Verify coherence between tcDMRG and FCI.
 * Here we take a conventional Hamiltonian, we encode it in the transcorrelated
 * format, and check that the energy is coherent with that of a "poor-man"
 * implementation of CISD.
 * Note that we do the check for both the conventional and the transcorrelated format.
 */
BOOST_FIXTURE_TEST_CASE(TestTCMolecular_H2_VersusFullCI, TranscorrelatedFixture)
{
    for (auto& iParameter: std::vector<DmrgParameters>{parametersH2Conventional_ConventionalFormat, parametersH2Conventional_TranscorrelatedFormat}) {
        iParameter.set("nsweeps", 10);
        iParameter.set("max_bond_dimension", 100);
        maquis::DMRGInterface<double> interface(iParameter);
        interface.optimize();
        auto energyDMRG = interface.energy();
        // Hand-made Full-CI
        auto lattice = Lattice(iParameter);
        auto model = Model<matrix, TwoU1>(lattice, iParameter);
        auto mpo = make_mpo(lattice, model);
        std::vector<MPS<matrix, TwoU1>> vectorOfMPS;
        for (const auto& iString: fullCIDeterminantsH2) {
            iParameter.set("init_state", "hf");
            iParameter.set("hf_occ", iString);
            vectorOfMPS.push_back(MPS<matrix, TwoU1>(lattice.size(), *(model.initializer(lattice, iParameter))));
        }
        matrix hamiltonianMatrix(vectorOfMPS.size(), vectorOfMPS.size(), 0.0);
        matrix eigenVectors(vectorOfMPS.size(), vectorOfMPS.size(), 0.0);
        alps::numeric::vector<double> eigenValues(vectorOfMPS.size(), 0.0);
        for (int iRow = 0; iRow < vectorOfMPS.size(); iRow++)
            for (int iCol = 0; iCol < vectorOfMPS.size(); iCol++)
                hamiltonianMatrix(iRow, iCol) = expval(vectorOfMPS[iRow], vectorOfMPS[iCol], mpo)/std::sqrt(norm(vectorOfMPS[iRow])*norm(vectorOfMPS[iCol]));
        alps::numeric::syev(hamiltonianMatrix, eigenVectors, eigenValues);
        BOOST_CHECK_CLOSE(eigenValues[vectorOfMPS.size()-1], energyDMRG, 1.0E-8);
    }
}

#endif // HAVE_TwoU1