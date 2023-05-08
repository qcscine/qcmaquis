/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/BenzeneFixture.h"
#include "dmrg/utils/time_stopper.h"
#include "dmrg/optimize/optimize.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
> symmetries;

/**
 * @brief Checks that [mpo_times_mps] gives results that are coherent with expval.
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE( Test_MPO_Times_MPS_ExpVal, S, symmetries, BenzeneFixture ) 
{
    // Generates the HF MPS
    auto lattice = Lattice(parametersBenzene);
    auto modelHF = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsHF = MPS<matrix, S>(lattice.size(), *(modelHF.initializer(lattice, parametersBenzene)));
    auto mpo = make_mpo(lattice, modelHF);
    // Calculates the MPS-MPO contraction
    auto traitClass = MPOTimesMPSTraitClass<matrix, S>(mpsHF, modelHF, lattice, 
                                                       modelHF.total_quantum_numbers(parametersBenzene),
                                                       parametersBenzene["max_bond_dimension"]);
    auto outputMPS = traitClass.applyMPO(mpo);
    // Calculates the energy in two ways
    auto energyFromMPSTimesMPO = overlap(mpsHF, outputMPS)/norm(mpsHF) + mpo.getCoreEnergy();
    auto energyFromExpVal = expval(mpsHF, mpo)/norm(mpsHF);
    BOOST_CHECK_CLOSE(energyFromMPSTimesMPO, energyFromExpVal, 1.E-10);
};

/**
 * @brief Checks consistency of the optimization of the benzene cation.
 * In this test we optimize the ground state of the benzene cation in two ways:
 * 1) With a standard optimization
 * 2) By ionizing the innermost orbital and, then, running a TS optimization.
 * The two strategies should give the same energy
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE( Test_MPO_Times_MPS_Ionization, S, symmetries, BenzeneFixture ) 
{
    // Types declaration
    using opt_base_t = optimizer_base<matrix, S, storage::disk>;
    int nSweeps = 50;
    // Conventional calculation
    parametersBenzene.set("hf_occ", "4,4,2,1,1,1");
    parametersBenzene.set("u1_total_charge1", 2);
    parametersBenzene.set("u1_total_charge2", 3);
    parametersBenzene.set("symmetry", "2u1pg");
    parametersBenzene.set("nsweeps", 10);
    maquis::DMRGInterface<double> interfaceCation(parametersBenzene);
    interfaceCation.optimize();
    auto cationicEnergy = interfaceCation.energy();
    // "By-hand" ionize-then-optimize
    parametersBenzene.set("hf_occ", "4,4,4,1,1,1");
    parametersBenzene.set("u1_total_charge1", 3);
    parametersBenzene.set("u1_total_charge2", 3);
    // Add noise to "move" the optimization away from the local energy minimum
    parametersBenzene.set("twosite_truncation", "heev_truncation");
    parametersBenzene.set("nsweeps", nSweeps);
    parametersBenzene.set("ngrowsweeps", 10);
    parametersBenzene.set("nmainsweeps", 10);
    parametersBenzene.set("alpha_initial", 1.0E-6);
    parametersBenzene.set("alpha_main", 1.0E-10);
    parametersBenzene.set("alpha_final", 0.);
    auto lattice = Lattice(parametersBenzene);
    auto model = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsHF = MPS<matrix, S>(lattice.size(), *(model.initializer(lattice, parametersBenzene)));
    auto mpo = make_mpo(lattice, model);
    // Creates the destructor operator
    auto traitClass = MPOTimesMPSTraitClass<matrix, S>(mpsHF, model, lattice, 
                                                       model.total_quantum_numbers(parametersBenzene),
                                                       parametersBenzene["max_bond_dimension"]);
    auto ionizedMPS = traitClass.ionizeMPS(0, generate_mpo::IonizedOrbital::Up);
    // "By hand" optimization
    auto stop_callback = time_stopper(static_cast<double>(parametersBenzene["run_seconds"]));
    std::shared_ptr<opt_base_t> optimizer;
    optimizer.reset( new ts_optimize<matrix, S, storage::disk>(ionizedMPS, mpo, parametersBenzene, stop_callback, lattice, 0) );
    for (int sweep=0; sweep < nSweeps; ++sweep)
      optimizer->sweep(sweep);
    auto energyByHand = expval(ionizedMPS, mpo)/norm(ionizedMPS);
    // This check could be made stricter, but with more sweeps
    BOOST_CHECK_CLOSE(cationicEnergy, energyByHand, 1.0E-7);
};
