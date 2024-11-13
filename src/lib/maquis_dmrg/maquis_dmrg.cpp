/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#include "maquis_dmrg.h"
#include <complex>
#include <utility>

#include "dmrg/sim/symmetry_factory.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/sim/interface_sim.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/models/measurements/chementropy.h"

namespace maquis {
#if defined(HAVE_SU2U1PG)
using SU2U1grp = SU2U1PG;
using TwoU1grp = TwoU1PG;
#elif defined(HAVE_SU2U1)
typedef SU2U1 SU2U1grp;
typedef TwoU1 TwoU1grp;
#endif

template <class ScalarType>
struct simulation_traits {
  using shared_ptr =
      std::shared_ptr<abstract_interface_sim<tmatrix<ScalarType>>>;
  template <class SymmGroup>
  struct F {
    using type = interface_sim<tmatrix<ScalarType>, SymmGroup>;
  };
};

template <typename ScalarType, Hamiltonian HamiltonianType>
struct DMRGInterface<ScalarType, HamiltonianType>::Impl {
  using sim_ptr = typename simulation_traits<ScalarType>::shared_ptr;
  sim_ptr sim;

  Impl(sim_ptr sim_) : sim(std::move(sim_)) {};
  ~Impl() = default;
};

template <typename ScalarType, Hamiltonian HamiltonianType>
DMRGInterface<ScalarType, HamiltonianType>::DMRGInterface(DmrgParameters& parms_
)
    : parms(parms_),
      impl_(new Impl(::dmrg::symmetry_factory<simulation_traits<ScalarType>>(
          parms_, parms_
      ))){};

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::run(const std::string& sim_type
) {
  try {
    impl_->sim->run(sim_type);
  } catch (std::exception& e) {
    maquis::cerr << "Exception thrown!" << std::endl;
    maquis::cerr << e.what() << std::endl;
    exit(1);
  }
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::optimize() {
  try {
    impl_->sim->run("optimize");
  } catch (std::exception& e) {
    maquis::cerr << "Exception thrown!" << std::endl;
    maquis::cerr << e.what() << std::endl;
    exit(1);
  }
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::evolve() {
  try {
    impl_->sim->run("evolve");
  } catch (std::exception& e) {
    maquis::cerr << "Exception thrown!" << std::endl;
    maquis::cerr << e.what() << std::endl;
    exit(1);
  }
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::run_measure() {
  try {
    impl_->sim->run_measure();
  } catch (std::exception& e) {
    maquis::cerr << "Exception thrown!" << std::endl;
    maquis::cerr << e.what() << std::endl;
    exit(1);
  }
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::runInversePowerIteration() {
  try {
    impl_->sim->run("ipi");
  } catch (std::exception& e) {
    maquis::cerr << "Exception thrown!" << std::endl;
    maquis::cerr << e.what() << std::endl;
    exit(1);
  }
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::runFEAST() {
  try {
    impl_->sim->run("feast");
  } catch (std::exception& e) {
    maquis::cerr << "Exception thrown!" << std::endl;
    maquis::cerr << e.what() << std::endl;
    throw;
    // exit(1);
  }
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::runTranscorrelated() {
  try {
    impl_->sim->run("transcorrelated");
  } catch (std::exception& e) {
    maquis::cerr << "Exception thrown!" << std::endl;
    maquis::cerr << e.what() << std::endl;
    throw;
  }
}

template <typename ScalarType, Hamiltonian HamiltonianType>
ScalarType DMRGInterface<ScalarType, HamiltonianType>::energy() {
  return impl_->sim->get_energy();
}

template <typename ScalarType, Hamiltonian HamiltonianType>
ScalarType DMRGInterface<ScalarType, HamiltonianType>::energyFEAST(int iState) {
  return impl_->sim->getFEASTEnergy(iState);
}

template <typename ScalarType, Hamiltonian HamiltonianType>
ScalarType DMRGInterface<ScalarType, HamiltonianType>::getCICoefficient(
    std::string determinantString
) {
  return impl_->sim->getCICoefficient(determinantString);
}

template <typename ScalarType, Hamiltonian HamiltonianType>
results_collector&
DMRGInterface<ScalarType, HamiltonianType>::get_iteration_results() {
  return impl_->sim->get_iteration_results();
}

template <typename ScalarType, Hamiltonian HamiltonianType>
int DMRGInterface<ScalarType, HamiltonianType>::get_last_sweep() {
  return impl_->sim->get_last_sweep();
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::measure() {
  measurements_ = impl_->sim->measure_out();
}

// #ifdef TRANSCORR_INTEGRALSC
template <class V, Hamiltonian HamiltonianType>
void DMRGInterface<V, HamiltonianType>::update_tc_integrals(
    const tc_integral_map& integrals
) {
  impl_->sim->update_tc_integrals(integrals);
}
// #endif

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::update_integrals(
    const integral_map<ScalarType>& integrals
) {
  impl_->sim->update_integrals(integrals);
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::update_integrals(
    std::string fileName
) {
  impl_->sim->update_integrals(fileName);
}
template <typename ScalarType, Hamiltonian HamiltonianType>
std::string DMRGInterface<ScalarType, HamiltonianType>::fiedler_order(
    int n_states, const std::vector<std::vector<int>>& hf_occupations,
    std::string checkpoint_name
) {
  results_map_type tmp_measurements = measurements_;
  DmrgParameters tmp_parms = parms;
  std::string ordering =
      impl_->sim->get_fiedler_order(n_states, hf_occupations, checkpoint_name);
  measurements_ = tmp_measurements;
  parms = tmp_parms;
  std::cout << "Fiedler order: " << ordering << "\n";
  return ordering;
}

template <typename ScalarType, Hamiltonian HamiltonianType>
const typename DMRGInterface<ScalarType, HamiltonianType>::results_map_type&
DMRGInterface<ScalarType, HamiltonianType>::measurements() {
  if (measurements_.empty()) {
    measure();
  }
  // This is probably not going to work if we call optimize() several times
  // TODO: handle also these cases!
  return measurements_;
};

template <typename ScalarType, Hamiltonian HamiltonianType>
const typename DMRGInterface<
    ScalarType, HamiltonianType>::meas_with_results_type&
DMRGInterface<ScalarType, HamiltonianType>::mutinf() {
  return measurements().at("mutinf");
}

// TODO: This does not work for 2U1/2U1PG symmetry because "oneptdm" measurement
// is not recognised by the model! Fix the model to recognise it!
template <typename ScalarType, Hamiltonian HamiltonianType>
const typename DMRGInterface<
    ScalarType, HamiltonianType>::meas_with_results_type&
DMRGInterface<ScalarType, HamiltonianType>::onerdm() {
  return measurements().at("oneptdm");
}

template <typename ScalarType, Hamiltonian HamiltonianType>
const typename DMRGInterface<
    ScalarType, HamiltonianType>::meas_with_results_type&
DMRGInterface<ScalarType, HamiltonianType>::onespdm() {
  return measurements().at("oneptspdm");
}

template <typename ScalarType, Hamiltonian HamiltonianType>
const typename DMRGInterface<
    ScalarType, HamiltonianType>::meas_with_results_type&
DMRGInterface<ScalarType, HamiltonianType>::twordm() {
  return measurements().at("twoptdm");
}

template <typename ScalarType, Hamiltonian HamiltonianType>
const typename DMRGInterface<
    ScalarType, HamiltonianType>::meas_with_results_type&
DMRGInterface<ScalarType, HamiltonianType>::threerdm() {
  parms.set("MEASURE[3rdm]", 1);  // required for 3-RDM measurement
  return measurements().at("threeptdm");
}

template <typename ScalarType, Hamiltonian HamiltonianType>
const typename DMRGInterface<
    ScalarType, HamiltonianType>::meas_with_results_type&
DMRGInterface<ScalarType, HamiltonianType>::fourrdm() {
  parms.set("MEASURE[4rdm]", 1);  // required for 4-RDM measurement
  return measurements().at("fourptdm");
}

template <typename ScalarType, Hamiltonian HamiltonianType>
const typename DMRGInterface<
    ScalarType, HamiltonianType>::meas_with_results_type&
DMRGInterface<ScalarType, HamiltonianType>::getMeasurement(std::string measName
) {
  if (measurements().find(measName) == measurements().end()) {
    throw std::runtime_error("Measurement not available!");
  }
  return measurements().at(measName);
}

#define measure_and_save_rdm(N)                     \
  BaseParameters meas_parms = parms.measurements(); \
  parms.erase_measurements();                       \
  parms.set("MEASURE[" #N "rdm]", 1);               \
  impl_->sim->run_measure();                        \
  parms.erase_measurements();                       \
  parms << meas_parms

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::measure_and_save_3rdm() {
  measure_and_save_rdm(3);
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::measure_and_save_4rdm() {
  // Clear all unnecessary measurements before running 4-RDM measurement
  // FIXME: clearing parms here has NO EFFECT on the measurements! This has to
  // be changed in another way! For now the measurements are modified in
  // maquis_cinterface.cpp, but it won't work if DMRGInterface is called
  // directly! Back up measurements
  measure_and_save_rdm(4);
}

#undef measure_and_save_rdm

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::measure_and_save_trans3rdm(
    const std::string& bra_name
) {
  BaseParameters meas_parms = parms.measurements();
  parms.erase_measurements();
  parms.set("MEASURE[trans3rdm]", bra_name);
  impl_->sim->run_measure();
  parms.erase_measurements();
  parms << meas_parms;
}

template <typename ScalarType, Hamiltonian HamiltonianType>
DMRGInterface<ScalarType, HamiltonianType>::~DMRGInterface() = default;

template <typename ScalarType, Hamiltonian HamiltonianType>
ScalarType DMRGInterface<ScalarType, HamiltonianType>::overlap(
    const std::string& aux_mps_name
) {
  return impl_->sim->get_overlap(aux_mps_name);
}

template <typename ScalarType, Hamiltonian HamiltonianType>
void DMRGInterface<ScalarType, HamiltonianType>::dump_parameters(
    const std::string& file
) {
  std::ofstream fs(file);
  fs << parms;
}

// Explicit template instantiation
template class DMRGInterface<double>;
template class DMRGInterface<std::complex<double>>;
// Check if these should be there
template class DMRGInterface<double, Hamiltonian::PreBO>;
template class DMRGInterface<double, Hamiltonian::VibrationalNMode>;
template class DMRGInterface<double, Hamiltonian::VibrationalCanonical>;
template class DMRGInterface<double, Hamiltonian::Vibronic>;
template class DMRGInterface<double, Hamiltonian::Excitonic>;
}  // namespace maquis
