/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef INTERFACE_SIM_H
#define INTERFACE_SIM_H

#include <sys/stat.h>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "dmrg/models/measurements/chementropy.h"
#include "dmrg/MetaSweepSimulations/FEASTLauncher.h"
#include "dmrg/SweepBasedAlgorithms/SweepOptimizationTypeTrait.h"
#include "dmrg/SweepBasedAlgorithms/SweepSimulationFactory.h"
#include "dmrg/evolve/TimeEvolutionSweep.h"
#include "dmrg/models/MolecularHamiltonians/measure_transform.hpp"
#include "dmrg/models/measurements.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/sim/abstract_sim.h"
#include "dmrg/sim/sim.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/utils/archive.h"
#include "dmrg/utils/checks.h"
#include "dmrg/utils/results_collector.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "integral_interface.h"
#include "utils/io.hpp"
#include "utils/timings.h"
#include "utils/traits.hpp"

namespace detail {
  template<typename T>
  static std::vector<int> sort_vector(std::vector<T>& vec_to_sort) {
    // prepare sorting
    std::vector<int> order(vec_to_sort.size());
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(), order.end(),
      [&vec_to_sort](size_t i1, size_t i2) {
        return maquis::real(vec_to_sort[i1]) < maquis::real(vec_to_sort[i2]);
      }
    );
    return order;
  }

  /**
   * @brief convert elements of a vector to a string.
   *
   * Basically convert elements of vector to a string, separatted by ','
   * @tparam K the type of the vector elements
   * @param v the vector to convert
   * @return s the string from the vector
   */
  template <class K>
  inline std::string vector_tostring(const std::vector<K>& v) {
    std::string s;
    for (int i = 0; i < v.size(); i++) {
      s += std::to_string(v[i]) + ((i < v.size()-1) ? "," : "") ;
    }
    return s;
  }

  /**
   * @brief convenience function to create checkpoint names for different states
   *
   * @param pname the name of the checkpoint file
   * @param state the nth state
   * @return the name of the checkpoint file
   */
  inline std::string checkpoint_name(const std::string & pname, int state) {
    return pname + ".checkpoint_state." + std::to_string(state) + ".h5";
  }

  /**
   * @brief Calculate the Laplacian
   *
   * @param mutI a matrix representation of the mutual information
   * @return laplacian a matrix representation laplacian
   */
  template<class Matrix>
  Matrix get_laplacian(const Matrix & mutI) {
    // Ported to C++ from fiedler.py
    Matrix laplacian(mutI.num_rows(), mutI.num_cols(), 0.0);
    for (int i = 0; i < mutI.num_rows(); i++) {
      for (int j = 0; j < mutI.num_cols(); j++) {
        laplacian(i,i) += mutI(i,j);
      }
    }
    laplacian -= mutI;
    return laplacian;
  }
}


// The sim class for interface-based DMRG runs and measurements
template <class Matrix, class SymmGroup>
class interface_sim : public sim<Matrix, SymmGroup>,
                      public abstract_interface_sim<Matrix> {
  // Types definition
  using base = sim<Matrix, SymmGroup>;
  using interface_base = abstract_interface_sim<Matrix>;
  using FEASTLauncherType = FEASTLauncher<Matrix, SymmGroup>;
  using measurements_type = typename base::measurements_type;
  using meas_with_results_type =
      typename interface_base::meas_with_results_type;
  using MPSType = MPS<Matrix, SymmGroup>;
  using ModelType = Model<Matrix, SymmGroup>;
  using results_map_type = typename interface_base::results_map_type;
  using FactoryType = SweepSimulationFactory<Matrix, SymmGroup, storage::disk>;
  using RealType = typename maquis::traits::real_type<Matrix>::type;
  using status_type = typename base::status_type;
  // Class inheritance from the sim object
  using base::all_measurements;
  using base::checkpoint_simulation;
  using base::chkpfolder;
  using base::init_site;
  using base::init_sweep;
  using base::lat;
  using base::model;
  using base::mpo;
  using base::mps;
  using base::parms;
  using base::results_archive_path;
  using base::rfile;
  using base::stop_callback;
  using base::sweep_measurements;

 public:
  /**
   * @brief Class constructor
   * Note that the base class is here the [sim] object.
   * @param parms_ parameter container
   */
  explicit interface_sim(DmrgParameters& parms_)
      : base(parms_), last_sweep_(init_sweep - 1) {}

  /** @brief Runs a DMRG-based optimization */
  void run(const std::string& simulationType) override {
    Timer timer(simulationType);
    timer.begin();
    if (simulationType == "optimize") {
      this->runAlternatingLeastSquares(
          "optimize", parms["nsweeps"].template as<int>(),
          parms["conv_thresh"].template as<double>(), model, parms
      );
    } else if (simulationType == "evolve") {
      this->runAlternatingLeastSquares(
          "evolve", parms["nsweeps"].template as<int>(),
          parms["conv_thresh"].template as<double>(), model, parms
      );
      // this->evolve();
    } else if (simulationType == "solve_linear_system") {
      this->runAlternatingLeastSquares(
          "linear_system", parms["nsweeps"].template as<int>(),
          parms["conv_thresh"].template as<double>(), model, parms
      );
    } else if (simulationType == "ipi") {
      this->runInversePowerIteration();
    } else if (simulationType == "feast") {
      this->runFEASTSimulation();
    } else if (simulationType == "transcorrelated") {
      this->runTranscorrelated();
    }
    timer.end();
  }

  /** @brief Runs a FEAST simulation */
  void runFEASTSimulation() {
#ifdef DMRG_FEAST
    try {
      feastMPSs_ =
          FEASTLauncherType::runFEASTSimulation(parms, model, lat, mpo);
    } catch (std::exception& e) {
      throw;
    }
    if (feastMPSs_) {
      for (int iState = 0; iState < feastMPSs_->size(); ++iState) {
        std::string filename = "FEAST_" + std::to_string(iState);
        checkpoint_simulation(feastMPSs_->operator[](iState), 1, -1, filename);
        dumpParameters(filename);
      }
    }
#else
    throw std::runtime_error(
        "Activate the BUILD_DMRG_FEAST Cmake flag for using DMRG[FEAST]"
    );
#endif
  }

  /** @brief Runs a IPI-based simulation */
  void runInversePowerIteration() {
    // Exctracts all relevant parameters
    double energyConvergenceThreshold = parms["ipi_sweep_energy_threshold"];
    double overlapConvergenceThreshold = parms["ipi_sweep_overlap_threshold"];
    int numberOfSweepsPerSystem = parms["nsweeps"];
    int numberOfOuterIterations = parms["ipi_iterations"];
    typename Matrix::value_type shift = parms["ipi_shift"];
    double convThreshOfLinSystem = parms["conv_thresh"];
    maquis::cout << " ===================================================== "
                 << std::endl;
    maquis::cout << "   STARTING DMRG[INVERSE POWER ITERATION] SIMULATION   "
                 << std::endl;
    maquis::cout << " ===================================================== "
                 << std::endl;
    maquis::cout << std::endl;
    maquis::cout << " IPI energy convergence threshold:   "
                 << energyConvergenceThreshold << std::endl;
    maquis::cout << " IPI overlap convergence threshold:  "
                 << overlapConvergenceThreshold << std::endl;
    maquis::cout << " Maximum number of IPI iterations:   "
                 << numberOfOuterIterations << std::endl;
    maquis::cout << " Number of sweeps per linear system: "
                 << numberOfSweepsPerSystem << std::endl;
    maquis::cout << " Shift parameter:                    " << shift
                 << std::endl;
    maquis::cout << std::endl;
    // Prepares data structure where to store results
    std::vector<RealType> energiesForIPIIteration;
    int nIpiIterations = 0;
    bool convergedOuter = false;
    double previousEnergy = this->get_energy();
    energiesForIPIIteration.push_back(previousEnergy);
    auto mpsBackup = this->mps;
    // IPI macroiteration
    while (!convergedOuter) {
      double nextEnergy, energyDifference;
      this->runAlternatingLeastSquares(
          "linear_system", numberOfSweepsPerSystem, convThreshOfLinSystem,
          model, parms
      );
      nIpiIterations += 1;
      nextEnergy = this->get_energy();
      energiesForIPIIteration.push_back(nextEnergy);
      energyDifference = std::abs(nextEnergy - previousEnergy);
      auto mpsOverlap = overlap(mpsBackup, this->mps) /
                        std::sqrt(norm(mpsBackup) * norm(this->mps));
      auto precision = std::cout.precision();
      maquis::cout << " === RESULTS FOR THE " << nIpiIterations
                   << "-th iteration ===" << std::endl;
      std::cout.precision(10);
      maquis::cout << " - Energy difference for iteration = " << nIpiIterations
                   << " = " << energyDifference << std::endl;
      maquis::cout << " - MPS overlap with solution at previous iteration = "
                   << std::abs(mpsOverlap) << std::endl;
      maquis::cout << std::endl;
      std::cout.precision(precision);
      // Checks convergence and, if not reached, starts a new IPI iteration
      if (nIpiIterations == numberOfOuterIterations ||
          energyDifference < energyConvergenceThreshold ||
          std::abs(1. - std::abs(mpsOverlap)) < overlapConvergenceThreshold) {
        std::string message =
            (nIpiIterations == numberOfOuterIterations)
                ? " --> MAXIMUM NUMBER OF IPI ITERATIONS REACHED"
                : " --> CONVERGENCE REACHED";
        maquis::cout << message << std::endl;
        convergedOuter = true;
      } else {
        maquis::cout << " --> CONVERGENCE NOT REACHED, STARTS NEW ITERATION"
                     << std::endl;
        previousEnergy = nextEnergy;
        mpsBackup = this->mps;
      }
    }
  }

  /**
   * @brief Generic ALS-based optimization.
   * This routine solve a problem (that can be either a propagation, the
   * solution of the linear system, or the energy minimization) with the ALS
   * algorithm. At the end, dumps the MPS, the energy, as well as the
   * measurements that were requested to be done at each microiteration.
   */
  void runAlternatingLeastSquares(
      std::string simulationType, int nSweeps, double energyThreshold,
      const ModelType& inputModel, DmrgParameters& inputParameters
  ) {
    bool verbose = inputParameters["verbose"] > 0;
    // Reads in input parameters
    int meas_each = parms["measure_each"];
    int chkp_each = parms["chkp_each"];
    // -- Optimizer initialization --
    if (parms["optimization"] == "singlesite") {
      // optimizer.reset( new ss_optimize<Matrix, SymmGroup, storage::disk>
      //                 (mps, mpo, parms, stop_callback, lat, init_site) );
      factory_ = std::make_unique<FactoryType>(
          simulationType, SweepOptimizationType::SingleSite, mps, mpo,
          inputParameters, inputModel, base::lat
      );
    } else if (parms["optimization"] == "twosite") {
      // optimizer.reset( new ts_optimize<Matrix, SymmGroup, storage::disk>
      //                 (mps, mpo, parms, stop_callback, lat, init_site) );
      factory_ = std::make_unique<FactoryType>(
          simulationType, SweepOptimizationType::TwoSite, mps, mpo,
          inputParameters, inputModel, base::lat
      );
    } else {
      throw std::runtime_error("Don't know this optimizer");
    }
    // Retrieve the measurements that should be always done.
    auto always_measurements = this->iteration_measurements(init_sweep);
    auto firstEnergy = this->get_energy();
    energies_.push_back(firstEnergy);

    if (verbose) {
      maquis::cout << "Initial energy is: " << std::setprecision(15)
                   << firstEnergy << std::endl;
    }
    // Run the sweep-based simulation.
    try {
      bool converged = false;
      for (int sweep = init_sweep; sweep < nSweeps; ++sweep) {
        factory_->runSingleSweep(sweep);
        energies_.push_back(this->get_energy());
        storage::disk::sync();
        if ((sweep + 1) % meas_each == 0 || (sweep + 1) == nSweeps) {
          dumpParametersAndIterResults(sweep);
          dumpEnergy(sweep);
          if (!rfile().empty() && always_measurements.size() > 0) {
            this->measure(
                this->results_archive_path(sweep) + "/results/",
                always_measurements
            );
            // stop simulation if an energy threshold has been specified
            int prev_sweep = sweep - meas_each;
            if (prev_sweep >= 0) {
              converged = checkEnergyConvergence(energyThreshold);
            }
          }
        }
        // Do not check convergence for propagation, since energy should be
        // conserved by definition
        if (!(simulationType == "evolve" && parms["imaginary_time"] == "no")) {
          converged = checkEnergyConvergence(energyThreshold);
        }
        last_sweep_ = sweep;
        /// write checkpoint
        bool stopped = stop_callback() || converged;
        if (stopped || (sweep + 1) % chkp_each == 0 || (sweep + 1) == nSweeps) {
          checkpoint_simulation(mps, sweep, -1);
        }
        if (stopped) {
          if (converged) {
            maquis::cout << "-- ALS CONVERGED --" << std::endl;
          }
          break;
        }
      }
      if (!converged) {
        maquis::cout << "-- ALS TERMINATED -- MAXIMUM ITERATIONS " << nSweeps
                     << " REACHED!" << std::endl;
      }
    } catch (const dmrg::time_limit& e) {
      maquis::cout << e.what() << " checkpointing partial result." << std::endl;
      checkpoint_simulation(mps, e.sweep(), e.site());
      dumpParametersAndIterResults(e.sweep());
      dumpEnergy(e.sweep());
    }
  }

  /**
   * @brief Runs a transcorrelated DMRG calculation.
   *
   * This simulation is composed by two steps:
   *
   *  1) a first TI-DMRG calculation associated with the non-transcorrelated
   * form of the target Hamiltonian. 2) a second imaginary-time TD-DMRG
   * calculation, which starts from the MPS optimized in the previous step.
   */
  void runTranscorrelated() {
#ifdef DMRG_TRANSCORRELATED
    // Extracts the relevant parameters
    int nSweepsTI = parms["transcorrelated_nsweeps_TI"];
    int nSweepsTC = parms["transcorrelated_nsweeps_TC"];
    double energyThreshold = parms["conv_thresh"];
    // Prints header
    maquis::cout << std::endl;
    maquis::cout << " ============================== " << std::endl;
    maquis::cout << "   STARTING tcDMRG SIMULATION = " << std::endl;
    maquis::cout << " ============================== " << std::endl;
    maquis::cout << " Preliminary TI-DMRG Sweeps: " << nSweepsTI << std::endl;
    maquis::cout << " tcDMRG Sweeps: " << nSweepsTC << std::endl;
    maquis::cout << " Energy Convergence Threshold: " << energyThreshold
                 << std::endl;
    maquis::cout << std::endl;
    // Preliminary TI calculation
    if (nSweepsTI > 0) {
      maquis::cout << " == STARTING THE TI-DMRG OPTIMIZATION == " << std::endl;
      maquis::cout << std::endl;
      auto conventionalParameterContainer = parms;
      auto conventionalModel = ModelType(lat, conventionalParameterContainer);
      mpo = make_mpo(base::lat, conventionalModel);
      this->runAlternatingLeastSquares(
          "optimize", nSweepsTI, energyThreshold, conventionalModel,
          conventionalParameterContainer
      );
    }
    // iTD-DMRG propagation
    init_sweep = nSweepsTI;
    if (nSweepsTC > 0) {
      maquis::cout << std::endl;
      maquis::cout << " == STARTING THE IMAGINARY-TIME TD-DMRG PROPAGATION  == "
                   << std::endl;
      maquis::cout << std::endl;
      // Generates the transcorrelated model
      auto transcorrelatedParametersContainer = parms;
      transcorrelatedParametersContainer.set(
          "transcorrelated_hamiltonian", "yes"
      );
      transcorrelatedParametersContainer.set("imaginary_time", "yes");
      auto transcorrelatedModel =
          ModelType(lat, transcorrelatedParametersContainer);
      mpo = make_mpo(base::lat, transcorrelatedModel);
      this->runAlternatingLeastSquares(
          "evolve", init_sweep + nSweepsTC, energyThreshold,
          transcorrelatedModel, transcorrelatedParametersContainer
      );
    }
#else
    throw std::runtime_error(
        "Activate the [BUILD_TRANSCORRELATED_DMRG] Cmake flag before running a "
        "tcDMRG calculation."
    );
#endif  // DMRG_TRANSCORRELATED
  }

  /** @brief Runs a propagation calculation */
  // AB For now it's mostly copy-pasted from optimize, should be rewritten in a
  // cleaner way.
  /*
  void evolve()
  {
#ifdef DMRG_TD
    // Types definition
    using EvolverType = TimeEvolutionSweep<Matrix, SymmGroup, storage::disk>;
    using SSEvolverType = SingleSiteTimeEvolution<Matrix, SymmGroup,
storage::disk>; using TSEvolverType = TwoSiteTimeEvolution<Matrix, SymmGroup,
storage::disk>;
    // Reads in input parameters
    int meas_each = parms["measure_each"];
    int chkp_each = parms["chkp_each"];
    // Optimizer initialization
    std::shared_ptr<EvolverType> evolver;
    if (parms["optimization"] == "singlesite")
        evolver = std::make_shared<SSEvolverType>(mps, mpo, parms,
stop_callback, init_site); else if(parms["optimization"] == "twosite") evolver =
std::make_shared<TSEvolverType>(mps, mpo, parms, stop_callback, init_site); else
        throw std::runtime_error("Evolution modality not recognized");
    measurements_type always_measurements =
this->iteration_measurements(init_sweep); int nSweeps = parms["nsweeps"]; try {
      for (int sweep=init_sweep; sweep < nSweeps; ++sweep) {
        evolver->evolve_sweep(sweep);
        storage::disk::sync();
        // Do measurements on the MPS after the time evolution
        if ((sweep + 1) % meas_each == 0 || (sweep + 1) == nSweeps) {
          // Measure energy
          auto energy = evolver->get_energy();
          iteration_results_ = evolver->iteration_results();
          // Measure observables specified in 'always_measure'
          if (!rfile().empty()) {
            measurements_type always_measure =
this->iteration_measurements(sweep); if (!parms["ALWAYS_MEASURE"].empty())
                this->measure(this->results_archive_path(sweep) + "/results/",
always_measure);
            // Write iteration results
            {
              storage::archive ar(rfile(), "w");
              ar[this->results_archive_path(sweep) + "/results"] <<
evolver->iteration_results(); ar[this->results_archive_path(sweep) +
"/results/Energy/mean/value"] << std::vector<double>(1, energy);

            }
          }
        }
        // Dump data in
        last_sweep_ = sweep;
        bool stopped = stop_callback();
        if (stopped || (sweep + 1) % chkp_each == 0 || (sweep + 1) ==
parms["nsweeps"]) checkpoint_simulation(mps, sweep, -1); if (stopped) break;
      }
    }
    catch (dmrg::time_limit const& e) {
      maquis::cout << e.what() << " checkpointing partial result." << std::endl;
      checkpoint_simulation(mps, e.sweep(), e.site());
      dumpParametersAndIterResults(e.sweep());
    }
#endif // DMRG_TD
  }
  */

  /** @brief Runs a measurement calculation */
  void run_measure() override {
    // if (this->get_last_sweep() < 0)
    //     throw std::runtime_error("Tried to measure before a sweep");
    this->measure("/spectrum/results/", all_measurements);
    // MPO creation
    MPO<Matrix, SymmGroup> mpoc = mpo;
    if (parms["use_compressed"]) {
      mpoc.compress(1e-12);
    }
    double energy;
    // Measures the energy
    if (parms["MEASURE[Energy]"]) {
      energy =
          maquis::real(expval(mps, mpoc)) / maquis::real(overlap(mps, mps));
      maquis::cout << "Energy: " << energy << std::endl;
      maquis::cout << "MPS norm: " << maquis::real(overlap(mps, mps))
                   << std::endl;
      if (!rfile().empty()) {
        storage::archive ar(rfile(), "w");
        ar["/spectrum/results/Energy/mean/value"]
            << std::vector<double>(1, energy);
      }
    }
    // Measures the energy variance
    if (parms["MEASURE[EnergyVariance]"]) {
      if (!parms["MEASURE[Energy]"]) {
        energy = maquis::real(expval(mps, mpoc));
      }
      auto traitClass = MPOTimesMPSTraitClass<Matrix, SymmGroup>(
          mps, model, base::lat, model.total_quantum_numbers(parms),
          parms["max_bond_dimension"]
      );
      auto outputMPS = traitClass.applyMPO(mpoc);
      auto energy2 = maquis::real(overlap(outputMPS, outputMPS) / norm(mps));
      maquis::cout << "Energy^2: " << energy2 << std::endl;
      maquis::cout << "Variance: " << energy2 - energy * energy << std::endl;
      if (!rfile().empty()) {
        storage::archive ar(rfile(), "w");
        ar["/spectrum/results/Energy^2/mean/value"]
            << std::vector<double>(1, energy2);
        ar["/spectrum/results/EnergyVariance/mean/value"]
            << std::vector<double>(1, energy2 - energy * energy);
      }
    }
#if defined(HAVE_TwoU1) || defined(HAVE_TwoU1PG)
    if (!rfile().empty()) {
      BaseParameters parms_meas;
      parms_meas = parms.twou1_measurements();
      if (!parms_meas.empty()) {
        measure_transform<Matrix, SymmGroup>(
        )(rfile(), "/spectrum/results", base::lat, mps, parms_meas);
      }
    } else {
      throw std::runtime_error(
          "Transformed measurements not implemented yet without checkpoints"
      );
    }
#endif
  }

  results_map_type measure_out() override {
    results_map_type ret;
    // Do not measure before a sweep
    if (this->get_last_sweep() < 0) {
      throw std::runtime_error("Tried to measure before a sweep");
    }
    // Run all measurements and fill the result map
    for (auto&& meas : all_measurements) {
      ret[meas.name()] =
          measure_and_save<Matrix, SymmGroup>(rfile(), "/spectrum/results", mps)
              .meas_out(meas);
    }
    // Measurements that require SU2U1->2U1 transformation
#if defined(HAVE_TwoU1) || defined(HAVE_TwoU1PG)
    BaseParameters parms_meas;
    parms_meas = parms.twou1_measurements();
    if (!parms_meas.empty()) {
      // Obtain a map with transformed measurements
      results_map_type transformed_meas =
          measure_transform<Matrix, SymmGroup>().meas_out(
              base::lat, mps, parms_meas, rfile(), "/spectrum/results"
          );
      // Merge transformed measurements with the remaining results
      ret.insert(transformed_meas.begin(), transformed_meas.end());
    }
#endif
    return ret;
  }

  /** @brief Gets the energy for the mps that is stored in the sim object */
  /** if it is a feastMPS, return the feast energy of the zeroth feast state*/
  RealType get_energy() override {
    if (!feastMPSs_) {
      return maquis::real(expval(mps, mpo) / overlap(mps, mps));
    } else {
      return getFEASTEnergy(0);
    }
  }

  /** @brief Gets the FEAST eigenstates - throws an exception if FEAST is not
   * run */
  auto getFEASTEigenstates() {
    if (!feastMPSs_) {
      throw std::runtime_error(
          "FEAST eigenstate requested before running a FEAST simulation"
      );
    } else {
      return feastMPSs_;
    }
  }

  /** @brief Gets the FEAST energies -- throws an exception if FEAST is not run
   */
  RealType getFEASTEnergy(int iState) const override {
    if (!feastMPSs_) {
      throw std::runtime_error(
          "FEAST energy requested before running a FEAST simulation"
      );
    } else if (iState >= feastMPSs_->size()) {
      std::string errorMessage =
          "FEAST energy requested for the" + std::to_string(iState) +
          "-th state, but only " + std::to_string(feastMPSs_->size()) +
          " states are available";
      throw std::runtime_error(errorMessage);
    } else {
      return maquis::real(
          expval(feastMPSs_->operator[](iState), mpo) /
          norm(feastMPSs_->operator[](iState))
      );
    }
  }

  /**
   * @brief Method to extract a CI coefficient associated to a given
   * determinant.
   *
   * Note that the format in which the input must be given is the same as for
   * the [hf_occ] (for electronic problems) or the [basis_state_generic]
   * initializer for the vibrational case.
   *
   * @param determinantString string associated with the target determinant
   * @return overlap value
   */
  typename Matrix::value_type getCICoefficient(std::string determinantString
  ) override {
    auto modifiedParameters = parms;
    std::string initState =
        (parms["MODEL"] == "quantum_chemistry") ? "hf" : "basis_state_generic";
    modifiedParameters.set("init_type", initState);
    if (parms["MODEL"] == "quantum_chemistry") {
      modifiedParameters.set("hf_occ", determinantString);
    } else {
      modifiedParameters.set("init_basis_state", determinantString);
    }
    auto mpsOverlap =
        MPSType(lat.size(), *(model.initializer(lat, modifiedParameters)));
    return overlap(mpsOverlap, mps) / std::sqrt(norm(mpsOverlap) * norm(mps));
  }

  /**
   * @brief Updates the integral and regenerates the data that depends on it.
   *
   * @Note This function erases existing entries of
   * 'integral_file' and 'integrals' from parameters
   * and sets 'integrals_binary' as source of integrals.
   *
   * @param integrals Integral map providing new integrals.
   */
  void update_integrals(
      const chem::integral_map<typename Matrix::value_type>& integrals
  ) override {
    // integrals are set later anyways
    // deleting old ones should be okay
    if (parms.is_set("integral_file")) {
      parms.erase("integral_file");
    }
    if (parms.is_set("integrals")) {
      parms.erase("integrals");
    }
    // if (parms.is_set("integral_file") || parms.is_set("integrals"))
    //     throw std::runtime_error("updating integrals in the interface not
    //     supported yet in the FCIDUMP format");
    parms.set("integrals_binary", chem::serialize(integrals));
    // construct new model and mpo with new integrals
    // hope this doesn't give any memory leaks
    model = Model<Matrix, SymmGroup>(lat, parms);
    mpo = make_mpo(lat, model);
    // check if MPS is still OK
    maquis::checks::right_end_check(mps, model.total_quantum_numbers(parms));
    all_measurements = model.measurements();
    all_measurements << overlap_measurements<Matrix, SymmGroup>(parms);
  }
    void update_tc_integrals(const chem::TranscorrMap<typename Matrix::value_type>& integrals)
    {

        // integrals are set later anyways
        // deleting old ones should be okay
        if (parms.is_set("integral_file")) {
          parms.erase("integral_file");
        }
        if(parms.is_set("integrals")){
          parms.erase("integrals");
        }
        // if (parms.is_set("integral_file") || parms.is_set("integrals"))
        //     throw std::runtime_error("updating integrals in the interface not supported yet in the FCIDUMP format");
        parms.set("integrals_binary", chem::serialize(integrals));

        //std::cout << " parms are set (and ints are updated) -> " << std::endl;
        //std::cout << parms << std::endl;

        // construct new model and mpo with new integrals
        // hope this doesn't give any memory leaks
        model = Model<Matrix, SymmGroup>(lat, parms);
        mpo = make_mpo(lat, model);

        // check if MPS is still OK
        maquis::checks::right_end_check(mps, model.total_quantum_numbers(parms));

        all_measurements = model.measurements();
        all_measurements << overlap_measurements<Matrix, SymmGroup>(parms);
  }

  results_collector& get_iteration_results() override {
    // If iteration_results is empty, we didn't perform the sweep yet, but
    // possibly loaded the MPS from a checkpoint so we need to load also
    // iteration results
    if (iteration_results_.empty()) {
      // If we are not loading from a checkpoint, last_sweep is set to -1
      // so we need to return an empty iteration_results vector
      if (get_last_sweep() < 0) return iteration_results_;

      // otherwise, we are restarting but there's something wrong with the
      // checkpoint
      if (!rfile().empty()) {
        try  // Load the iteration results from the last sweep
        {
          storage::archive ar(rfile(), "r");
          ar[results_archive_path(last_sweep_) + "/results"] >>
              iteration_results_;
        } catch (std::exception& e) {
          maquis::cerr << e.what() << std::endl;
          throw std::runtime_error(
              "Error reading iteration results from checkpoint."
          );
        }
      } else
        throw std::runtime_error(
            "No result file specified for restart -- cannot read iteration "
            "results!"
        );
    }

    return iteration_results_;
  }
  /** @brief Updates the integrals from a new file */
  void update_integrals(std::string fileName) {
    parms.set("integral_file", fileName);
    model = Model<Matrix, SymmGroup>(lat, parms);
    mpo = make_mpo(lat, model);
    // check if MPS is still OK
    maquis::checks::right_end_check(mps, model.total_quantum_numbers(parms));
    all_measurements = model.measurements();
    all_measurements << overlap_measurements<Matrix, SymmGroup>(parms);
  }

  /** @brief Get the overlap of the MPS with another MPS, which is loaded from a
   * chkp file */
   typename Matrix::value_type get_overlap(
      const std::string& aux_filename
  ) override {
    maquis::checks::symmetry_check(parms, aux_filename);
    MPS<Matrix, SymmGroup> aux_mps;
    load(aux_filename, aux_mps);
    return overlap(aux_mps, this->mps);
  }

  /**
   * @brief Generate Fiedler ordering
   *
   * @param n_states int number of states
   * @param hf_occupations vector of vectors of ints the occupation for each state
   * @return order a string with the orbital order based on fiedler ordering
   */

  // template<typename std::enable_if<std::is_same<ScalarType, double>::value, int>::type = 0>
  // template<typename = typename std::enable_if <std::is_same<ScalarType, double>::value>::type>
  std::string get_fiedler_order(int n_states, const std::vector<std::vector<int>>& hf_occupations, std::string checkpoint_name) {
    maquis::cout << "-----------------------------------------------------------------" << std::endl;
    maquis::cout << "Start Fiedeler Ordering" << std::endl;
    maquis::cout << "-----------------------------------------------------------------" << std::endl;

    using ScalarType = typename Matrix::value_type;
    using meas_with_results_type = std::pair<std::vector<std::vector<int> >, std::vector<ScalarType> >;
    using results_map_type = std::map<std::string, meas_with_results_type>;
    // start new measurements
    parms.erase_measurements();

    if (!hf_occupations.empty()) {
      assert(hf_occupations.size() == n_states);
    }

    if (parms.is_set("orbital_order")) {
      // reset to default orbital order if some order is present
      // if this isn't done, there're strange side-effects
      int L = parms["L"];
      std::vector<int> v(L);
      std::iota(v.begin(), v.end(), 1);
      parms.set("orbital_order", detail::vector_tostring(v));
    }

    // we need mutual information for the Fiedler ordering
    parms.set("MEASURE[ChemEntropy]", 1);

    std::vector<results_map_type> measurements;
    measurements.reserve(n_states);

    // set sweeps and m, same values as in the old python interface
    parms.set("nsweeps", 4);
    // if (parms.is_set("init_bond_dimension")) {
    //   int init_bond_dimension = parms["init_bond_dimension"];
    //   parms.set("max_bond_dimension", init_bond_dimension);
    // } else {
    if (parms.is_set("L")) {
      parms.set("max_bond_dimension", parms["L"] > 24 ? 256 : 128);
    } else {
      throw std::runtime_error("L not defined for a starting guess calculation!");
    }
    // }

    //if(parms.is_set("feast_num_states")) {
    //    measurements.reserve(parms["feast_num_states"]);
    //    this->run("feast");
    //    // auto feast_mps = this->getFEASTEigenstates();
    //    for (int iState = 0; iState < feastMPSs_->size(); ++iState) {
    //      auto& mps = (*feastMPSs_)[iState];
    //      results_map_type ret;
    //      for (auto&& meas: all_measurements){
    //        ret[meas.name()] = measure_and_save<Matrix,SymmGroup>(rfile(), "/spectrum/results", mps).meas_out(meas);
    //        measurements.emplace_back(std::move(this->measure_out()));
    //      }
    //    }

    //} else {
    // Do it for each state
    for (int i = 0; i < n_states; i++) {
      // set correct checkpoints and result file names
      std::string chkpfile = detail::checkpoint_name(checkpoint_name, i);
      parms.set("chkpfile", chkpfile);

      // set HF occupation
      if (!hf_occupations.empty()) {
        parms.set("hf_occ", detail::vector_tostring(hf_occupations[i]));
      }

      // if excited state
      if (i > 0) {
        maquis::cout << "excited states" << std::endl;
        parms.set("n_ortho_states", i-1);
        std::string all_ortho_states;
        for (int j = 0; j < i; j++) {
          maquis::cout << "add state: " << j << std::endl;
          all_ortho_states += detail::checkpoint_name(checkpoint_name, j) + ((j < i-1) ? " " : "");
        }
        std::cout << all_ortho_states << std::endl;
        parms.set("ortho_states", all_ortho_states);
      }

      // do dmrg calculation
      maquis::cout << "Optimize for Fiedler" << std::endl;
      this->run("optimize");
      measurements.emplace_back(std::move(this->measure_out()));
    }
    // }

    // Get state-average single-orbital entropies and mutual information
    // Collect mutual information from all the states
    std::vector<Matrix> mutI;

    // Calculate S1 only if CI-DEAS is requested and mutual information only if Fiedler ordering is requested
    std::vector<Matrix> s1_;
    s1_.reserve(n_states);
    mutI.reserve(n_states);

    for (int i = 0; i < n_states; i++) {
      // get the entropy data
      EntanglementData<Matrix> em(measurements[i]);

      s1_.emplace_back(std::move(em.s1()));
      mutI.emplace_back(std::move(em.I()));
    }

    // Calculate average mutual information
    Matrix SAmutI(mutI[0].num_rows(), mutI[0].num_cols(),0.0);
    for (auto& n : mutI) {
      SAmutI += n;
    }

    // Divide mutual information by the number of states: irrelevant for Fiedler ordering
    // but let's still do it for the consistency
    //SAmutI /= n_states;
    Matrix SA_mutI_ = SAmutI;


    // TODO: implement also Block fiedler ordering per symmetry

    // get Laplacian of the average mutual information
    Matrix L = detail::get_laplacian(SA_mutI_);

    if (L.num_rows() < 2) {
      throw std::runtime_error("Fiedler vector orbital ordering doesn't work for only one orbital!");
    }

    // get eigenvectors and eigenvalues of the Laplacian
    Matrix evecs(L.num_rows(), L.num_cols());
    std::vector<double> evals(L.num_rows());
    alps::numeric::syev(L,evecs,evals);

    // get the Fiedler vector, i.e. the eigenvector corresponding to the second lowest eigenvalue of the Laplacian
    // The eigenvalues in evecs are assumed to be sorted starting from the highest eigenvalue
    // i.e. the second lowest eigenvalue has an index L-2
    auto fv_col = evecs.col(L.num_rows()-2);
    std::vector<ScalarType> fiedler_vector(fv_col.first, fv_col.second);

    // old
    /*
      // prepare ordering. first create a vector with indices 0..L-1 in ascending order
      std::vector<int> order(fiedler_vector.size());
      std::iota(order.begin(), order.end(), 0);

      // Sort the order vector according to the Fiedler vector
      std::sort(order.begin(), order.end(),
          [&fiedler_vector](size_t i1, size_t i2) {
            return fiedler_vector[i1] < fiedler_vector[i2];
          }
      );
    */
    std::vector<int> order  = detail::sort_vector(fiedler_vector);

    // add 1 to each element because in the parameters our counting starts with 1
    // This is not used right?
    // for (auto&& n: order) n++;
    // std::transform(order.begin(), order.end(), order.begin(), [](int i){ return i+1; });

    // convert the ordering into a string
    maquis::cout << "-----------------------------------------------------------------" << std::endl;
    maquis::cout << "End Fiedeler Ordering" << std::endl;
    maquis::cout << "-----------------------------------------------------------------" << std::endl;
    return detail::vector_tostring(order);
  }

  /** @brief Getter for the number of sweeps that have been run */
  int get_last_sweep() override { return last_sweep_; }

  /** @brief Class destructor */
  ~interface_sim() { storage::disk::sync(); }

 private:
  /** @brief Dumps the parameters and the iteration results to a result file */
  void dumpParametersAndIterResults(int iSweep) {
    // iteration_results_ = optimizer->iteration_results();
    iteration_results_ = factory_->getIterationResults();
    /// write iteration results if result files are specified
    if (!rfile().empty()) {
      storage::archive ar(rfile(), "w");
      ar[results_archive_path(iSweep) + "/parameters"] << parms;
      ar[results_archive_path(iSweep) + "/results"] << iteration_results_;
    }
  }

  /** @brief Dumps the energy to the result file */
  void dumpEnergy(int iSweep) {
    auto energy = this->get_energy();
    energies_.push_back(energy);
    if (!rfile().empty()) {
      auto energy = this->get_energy();
      storage::archive ar(rfile(), "w");
      ar[this->results_archive_path(iSweep) + "/results/Energy/mean/value"]
          << std::vector<double>(1, energy);
    }
  }

  /**  @brief Checks energy convergence of the sweep-based optimization */
  bool checkEnergyConvergence(double convergenceThreshold) {
    if (energies_.size() < 2) {  // Not yet sufficient number of iterations
      return false;
    }
    auto eDiff = std::abs(*(energies_.end() - 2) - *(energies_.end() - 1));
    return (eDiff < convergenceThreshold);
  }

  /** @brief Returns the path where the result of a given sweep are stored */
  std::string results_archive_path(int sweep) const {
    status_type status;
    status["sweep"] = sweep;
    return base::results_archive_path(status);
  }

  /** @brief Dumps the simulation results to the checkpoint file */
  void checkpoint_simulation(
      const MPS<Matrix, SymmGroup>& state, int sweep, int site,
      std::string filename = ""
  ) {
    status_type status;
    status["sweep"] = sweep;
    status["site"] = site;
    return base::checkpoint_simulation(state, status, filename);
  }

  void dumpParameters(std::string filename = "") {
    if (!chkpfolder().empty()) {
      std::string chkpfilename;
      if (filename.empty()) {
        chkpfilename = chkpfolder();
      } else {
        chkpfilename = chkpfolder() + "_" + filename;
      }
      storage::archive ar(chkpfilename + "/props.h5", "w");
      ar["/parameters"] << parms;
    }
  }

  // +-- Class members --+
  results_collector iteration_results_;
  int last_sweep_;
  std::unique_ptr<FactoryType> factory_;
  std::vector<RealType> energies_;
  std::shared_ptr<std::vector<MPSType>> feastMPSs_;
};

#endif
