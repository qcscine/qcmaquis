/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021- Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MODULE MPSOverlapElectronic

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/BenzeneFixture.h"
#include "Fixtures/RelativisticFixture.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** @brief Checks that overlap and norm give coherent results for the electronic case */
BOOST_FIXTURE_TEST_CASE_TEMPLATE( Test_MPS_Overlap_Electronic, S, symmetries, BenzeneFixture )
{
    // Modifies the last parameters
    auto lattice = Lattice(parametersBenzene);
    auto model = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsHF = MPS<matrix, S>(lattice.size(), *(model.initializer(lattice, parametersBenzene)));
    double normHF = norm(mpsHF);
    double overlapHF = overlap(mpsHF, mpsHF);
    BOOST_CHECK_CLOSE(normHF, overlapHF, 1.E-10);
}

/** @brief Checks Hermitianity of the overlap calculation for the real-valued electronic case */
BOOST_FIXTURE_TEST_CASE_TEMPLATE( Test_MPS_Overlap_Hermitian_Electronic, S, symmetries, BenzeneFixture )
{
    // Global variables
    auto lattice = Lattice(parametersBenzene);
    auto model = Model<matrix, S>(lattice, parametersBenzene);
    // Modifies the init parameters to create two different MPSs
    parametersBenzene.set("init_state", "const");
    auto mpsConst = MPS<matrix, S>(lattice.size(), *(model.initializer(lattice, parametersBenzene)));
    parametersBenzene.set("init_state", "default");
    auto mpsDefault = MPS<matrix, S>(lattice.size(), *(model.initializer(lattice, parametersBenzene)));
    // Calculates the overlap
    double overlapOriginal = overlap(mpsConst, mpsDefault);
    double overlapHerm = overlap(mpsDefault, mpsConst);
    BOOST_CHECK_CLOSE(overlapOriginal, overlapHerm, 1.E-10);
}

#ifdef HAVE_U1DG

/** @brief Checks Hermitianity of the overlap calculation for the relativistic (complex-valued) electronic case */
BOOST_FIXTURE_TEST_CASE( Test_MPS_Overlap_Hermitian_Electronic_Complex, RelativisticFixture )
{
    // Global variables
    auto lattice = Lattice(parametersComplex);
    auto model = Model<cmatrix, U1DG>(lattice, parametersComplex);
    // Modifies the init parameters to create two different MPSs
    parametersComplex.set("init_state", "const");
    auto mpsConst = MPS<cmatrix, U1DG>(lattice.size(), *(model.initializer(lattice, parametersComplex)));
    parametersComplex.set("init_state", "default");
    auto mpsDefault = MPS<cmatrix, U1DG>(lattice.size(), *(model.initializer(lattice, parametersComplex)));
    // Calculates the overlap
    std::complex<double> overlapOriginal = overlap(mpsConst, mpsDefault);
    std::complex<double> overlapHerm = overlap(mpsDefault, mpsConst);
    BOOST_CHECK_CLOSE(std::real(overlapOriginal), std::real(overlapHerm), 1.E-10);
    BOOST_CHECK_CLOSE(std::imag(overlapOriginal), -std::imag(overlapHerm), 1.E-10);
}

#endif
