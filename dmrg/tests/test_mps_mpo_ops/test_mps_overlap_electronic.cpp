/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
    parametersBenzene.set("init_type", "const");
    auto mpsConst = MPS<matrix, S>(lattice.size(), *(model.initializer(lattice, parametersBenzene)));
    parametersBenzene.set("init_type", "default");
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
    parametersComplex.set("init_type", "const");
    auto mpsConst = MPS<cmatrix, U1DG>(lattice.size(), *(model.initializer(lattice, parametersComplex)));
    parametersComplex.set("init_type", "default");
    auto mpsDefault = MPS<cmatrix, U1DG>(lattice.size(), *(model.initializer(lattice, parametersComplex)));
    // Calculates the overlap
    std::complex<double> overlapOriginal = overlap(mpsConst, mpsDefault);
    std::complex<double> overlapHerm = overlap(mpsDefault, mpsConst);
    BOOST_CHECK_CLOSE(std::real(overlapOriginal), std::real(overlapHerm), 1.E-10);
    BOOST_CHECK_CLOSE(std::imag(overlapOriginal), -std::imag(overlapHerm), 1.E-10);
}

#endif
