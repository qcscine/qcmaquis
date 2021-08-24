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

#define BOOST_TEST_MODULE TimeEvolutionSweep

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "utils/fpcomparison.h"
#include <iostream>
#include "dmrg/utils/time_stopper.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/evolve/TimeEvolutionSweep.h"
#include "Fixtures/TimeEvolversFixture.h"

// Note that the relativistic groups are managed separately
typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** Checks the constructor of a site shifter object */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestConstructorSingleSiteEvolution, S, symmetries, TestTimeEvolverFixture) {
    using SingleSiteTimeEvolution = SingleSiteTimeEvolution<cmatrix, S, storage::disk>;
    // In both cases, we construct everything as is done in [sim.hpp]
    auto lat = Lattice(parametersH2FourOrbitals);
    auto model = Model<cmatrix, S>(lat, parametersH2FourOrbitals);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, parametersH2FourOrbitals)));
    auto latticeSize = mpo.length();
    mps.normalize_right();
    // Prepares the boundaries
    time_stopper stop_callback(10000.);
    auto ssEvolver = SingleSiteTimeEvolution(mps, mpo, parametersH2FourOrbitals, stop_callback);
    auto initialEnergy = expval(mps, mpo);
    for (int iSweep = 0; iSweep < 10; iSweep++)
        ssEvolver.evolve_sweep(iSweep);
    auto finalEnergy = expval(mps, mpo);
    BOOST_CHECK_CLOSE(maquis::real(initialEnergy), maquis::real(finalEnergy), 1.0E-10);
}

/** Checks the constructor of a site shifter object */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestConstructorTwoSiteEvolution, S, symmetries, TestTimeEvolverFixture) {
    using TwoSiteTimeEvolution = TwoSiteTimeEvolution<cmatrix, S, storage::disk>;
    // We construct everything as is done in [sim.hpp]
    auto lat = Lattice(parametersH2FourOrbitals);
    auto model = Model<cmatrix, S>(lat, parametersH2FourOrbitals);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, parametersH2FourOrbitals)));
    auto latticeSize = mpo.length();
    mps.normalize_right();
    // Prepares the boundaries
    time_stopper stop_callback(10000.);
    auto tsEvolver = TwoSiteTimeEvolution(mps, mpo, parametersH2FourOrbitals, stop_callback);
    auto initialEnergy = expval(mps, mpo);
    for (int iSweep = 0; iSweep < 10; iSweep++)
        tsEvolver.evolve_sweep(iSweep);
    auto finalEnergy = expval(mps, mpo);
    BOOST_CHECK_CLOSE(maquis::real(initialEnergy), maquis::real(finalEnergy), 1.0E-10);
}


BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestSingleSiteImaginaryTimeVsTI, S, symmetries, TestTimeEvolverFixture) {
    using SingleSiteTimeEvolution = SingleSiteTimeEvolution<cmatrix, S, storage::disk>;
    // We construct everything as is done in [sim.hpp]
    auto lat = Lattice(parametersH2FourOrbitalsImaginary);
    auto model = Model<cmatrix, S>(lat, parametersH2FourOrbitalsImaginary);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, parametersH2FourOrbitalsImaginary)));
    auto latticeSize = mpo.length();
    mps.normalize_right();
    // Prepares the boundaries
    time_stopper stop_callback(10000.);
    auto ssEvolver = SingleSiteTimeEvolution(mps, mpo, parametersH2FourOrbitalsImaginary, stop_callback);
    // We extract the number of sweeps from the parameters to ensure consistency with TI-DMRG
    int nSweeps = parametersH2FourOrbitalsImaginary["nsweeps"];
    for (int iSweep = 0; iSweep < nSweeps; iSweep++)
        ssEvolver.evolve_sweep(iSweep);
    auto iTDEnergy = expval(mps, mpo);
    // Calculates the same energy via conventional TI-DMRG
    auto parameterCopy = parametersH2FourOrbitalsReal;
    parameterCopy.set("symmetry", returnStringRepresentation<S>());
    maquis::DMRGInterface<double> interface(parameterCopy);
    interface.optimize();
    // Checks consistency
    BOOST_CHECK_CLOSE(maquis::real(iTDEnergy), interface.energy(), 1.0E-10);
}


BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestTwoSiteImaginaryTimeVsTI, S, symmetries, TestTimeEvolverFixture) {
    using TwoSiteTimeEvolution = TwoSiteTimeEvolution<cmatrix, S, storage::disk>;
    // We construct everything as is done in [sim.hpp]
    auto lat = Lattice(parametersH2FourOrbitalsImaginary);
    auto model = Model<cmatrix, S>(lat, parametersH2FourOrbitalsImaginary);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, parametersH2FourOrbitalsImaginary)));
    auto latticeSize = mpo.length();
    mps.normalize_right();
    // Prepares the boundaries
    time_stopper stop_callback(10000.);
    auto tsEvolver = TwoSiteTimeEvolution(mps, mpo, parametersH2FourOrbitalsImaginary, stop_callback);
    // We extract the number of sweeps from the parameters to ensure consistency with TI-DMRG
    int nSweeps = parametersH2FourOrbitalsImaginary["nsweeps"];
    for (int iSweep = 0; iSweep < nSweeps; iSweep++)
        tsEvolver.evolve_sweep(iSweep);
    auto iTDEnergy = expval(mps, mpo);
    // Calculates the same energy via conventional TI-DMRG
    auto parameterCopy = parametersH2FourOrbitalsReal;
    parameterCopy.set("symmetry", returnStringRepresentation<S>());
    maquis::DMRGInterface<double> interface(parameterCopy);
    interface.optimize();
    // Checks consistency
    BOOST_CHECK_CLOSE(maquis::real(iTDEnergy), interface.energy(), 1.0E-10);
}

#ifdef HAVE_U1DG

/** 
 * @brief Real-time SS for relativistic Hamiltonians
 * Here we check that after 10 sweeps the energy is conserved.
 */
BOOST_FIXTURE_TEST_CASE(TestTwoSiteRealTimeVsTIRelativistic, TestTimeEvolverFixture) {
    using S = U1DG;
    using TwoSiteTimeEvolution = TwoSiteTimeEvolution<cmatrix, S, storage::disk>;
    // We construct everything as is done in [sim.hpp]
    auto lat = Lattice(parametersRelativistic);
    auto model = Model<cmatrix, S>(lat, parametersRelativistic);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, parametersRelativistic)));
    auto latticeSize = mpo.length();
    mps.normalize_right();
    // We modify the parameters taken from the fixture class, which are set on the 
    // imaginary-time case.
    parametersRelativistic.set("imaginary_time", "no");
    parametersRelativistic.set("time_units", "as");
    parametersRelativistic.set("TD_backpropagation", "yes");
    // Prepares the boundaries
    time_stopper stop_callback(10000.);
    auto tsEvolver = TwoSiteTimeEvolution(mps, mpo, parametersRelativistic, stop_callback);
    auto energyBefore = expval(mps, mpo)/norm(mps);
    for (int iSweep = 0; iSweep < 20; iSweep++)
        tsEvolver.evolve_sweep(iSweep);
    auto energyAfter = expval(mps, mpo)/norm(mps);
    // Checks consistency
    BOOST_CHECK_CLOSE(maquis::real(energyBefore), maquis::real(energyAfter), 1.0E-8);
}

/** 
 * @brief Imaginary-time TS for relativistic Hamiltonians 
 * Note that the reference energy is taken from [test_rel]
 */
BOOST_FIXTURE_TEST_CASE(TestSingleSiteImaginaryTimeVsTIRelativistic, TestTimeEvolverFixture) {
    using S = U1DG;
    using TwoSiteTimeEvolution = TwoSiteTimeEvolution<cmatrix, S, storage::disk>;
    // We construct everything as is done in [sim.hpp]
    auto lat = Lattice(parametersRelativistic);
    auto model = Model<cmatrix, S>(lat, parametersRelativistic);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, parametersRelativistic)));
    auto latticeSize = mpo.length();
    mps.normalize_right();
    // Prepares the boundaries
    time_stopper stop_callback(10000.);
    auto tsEvolver = TwoSiteTimeEvolution(mps, mpo, parametersRelativistic, stop_callback);
    // We extract the number of sweeps from the parameters to ensure consistency with TI-DMRG
    int nSweeps = parametersRelativistic["nsweeps"];
    for (int iSweep = 0; iSweep < nSweeps; iSweep++)
        tsEvolver.evolve_sweep(iSweep);
    auto iTDEnergy = expval(mps, mpo);
    // Checks consistency
    BOOST_CHECK_CLOSE(maquis::real(iTDEnergy), -1.0780470133e+02, 1.0E-7);
}

#endif