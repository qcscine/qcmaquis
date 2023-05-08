/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef INTERFACE_SIM_H
#define INTERFACE_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "dmrg/sim/sim.h"
// #include "dmrg/optimize/optimize.h"
#include "dmrg/evolve/TimeEvolutionSweep.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/models/chem/measure_transform.hpp"
#include "integral_interface.h"
#include "dmrg/utils/results_collector.h"
#include "dmrg/MetaSweepSimulations/FEASTLauncher.h"
#include "dmrg/SweepBasedAlgorithms/SweepSimulationFactory.h"

// The sim class for interface-based DMRG runs and measurements
template <class Matrix, class SymmGroup>
class interface_sim : public sim<Matrix, SymmGroup>, public abstract_interface_sim<Matrix> {
  // Types definition
  using base = sim<Matrix, SymmGroup>;
  using interface_base = abstract_interface_sim<Matrix>;
  using FEASTLauncherType = FEASTLauncher<Matrix, SymmGroup>;
  using measurements_type = typename base::measurements_type;
  using meas_with_results_type = typename interface_base::meas_with_results_type;
  using MPSType = MPS<Matrix, SymmGroup>;
  using results_map_type = typename interface_base::results_map_type;
  using FactoryType = SweepSimulationFactory<Matrix, SymmGroup, storage::disk>;
  using RealType = typename maquis::traits::real_type<Matrix>::type;
  using status_type = typename base::status_type;
  // Class inheritance from the sim object
  using base::mps;
  using base::mpo;
  using base::parms;
  using base::all_measurements;
  using base::sweep_measurements;
  using base::stop_callback;
  using base::init_sweep;
  using base::init_site;
  using base::rfile;
  using base::chkpfolder;
  using base::lat;
  using base::model;

public:
  /**
   * @brief Class constructor
   * Note that the base class is here the [sim] object.
   * @param parms_ parameter container
   */
  explicit interface_sim(DmrgParameters & parms_) : base(parms_), last_sweep_(init_sweep-1) { }

  /** @brief Runs a DMRG-based optimization */
  void run(const std::string& simulationType) {
    if (simulationType == "optimize")
      this->runAlternatingLeastSquares("optimize", parms["nsweeps"].template as<int>(), parms["conv_thresh"].template as<double>());
    else if (simulationType == "evolve")
      this->runAlternatingLeastSquares("evolve", parms["nsweeps"].template as<int>(), parms["conv_thresh"].template as<double>());
      //this->evolve();
    else if (simulationType == "solve_linear_system")
      this->runAlternatingLeastSquares("linear_system", parms["nsweeps"].template as<int>(), parms["conv_thresh"].template as<double>());
    else if (simulationType == "ipi")
      this->runInversePowerIteration();
    else if (simulationType == "feast")
      this->runFEASTSimulation();
  }

  /** @brief Runs a FEAST simulation */
  // TODO: fix the MPS that is actually extracted -- it should be not necessarily th 0-th one.
  void runFEASTSimulation() {
#ifdef DMRG_FEAST
    try {
      feastMPSs_ = FEASTLauncherType::runFEASTSimulation(parms, model, lat, mpo);
    }
    catch (std::exception& e) {
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
    throw std::runtime_error("Activate the BUILD_DMRG_FEAST Cmake flag for using DMRG[FEAST]");
#endif
  }

  /** @brief Runs a IPI-based simulation */
  void runInversePowerIteration() {
    // Exctracts all relevant parameters
    double energyConvergenceThreshold = parms["ipi_sweep_energy_threshold"];
    double overlapConvergenceThreshold = parms["ipi_sweep_overlap_threshold"];
    int numberOfSweepPerSystem = parms["ipi_sweeps_per_system"];
    int numberOfOuterIterations = parms["ipi_iterations"];
    typename Matrix::value_type shift = parms["ipi_shift"];
    maquis::cout << " ===================================================== " << std::endl;
    maquis::cout << "   STARTING DMRG[INVERSE POWER ITERATION] SIMULATION = " << std::endl;
    maquis::cout << " ===================================================== " << std::endl;
    maquis::cout << std::endl;
    maquis::cout << " IPI energy convergence threshold:   " << energyConvergenceThreshold << std::endl;
    maquis::cout << " IPI overlap convergence threshold:  " << overlapConvergenceThreshold << std::endl;
    maquis::cout << " Number of sweeps per linear system: " << numberOfSweepPerSystem << std::endl;
    maquis::cout << " Shift parameter: " << shift << std::endl;
    maquis::cout << std::endl;
    // Prepares data structure where to store results
    std::vector<RealType> energiesForIPIIteration;
    int nIpiIterations=0;
    bool convergedOuter=false;
    double previousEnergy = this->get_energy();
    energiesForIPIIteration.push_back(previousEnergy);
    auto mpsBackup = this->mps;
    // IPI macroiteration
    while (!convergedOuter) {
      double nextEnergy, energyDifference;
      this->runAlternatingLeastSquares("linear_system", numberOfSweepPerSystem, 0.);
      nIpiIterations += 1;
      nextEnergy = this->get_energy();
      energiesForIPIIteration.push_back(nextEnergy);
      energyDifference = std::fabs(nextEnergy - previousEnergy);
      auto mpsOverlap = overlap(mpsBackup, this->mps)/std::sqrt(norm(mpsBackup)*norm(this->mps));
      auto precision = std::cout.precision();
      maquis::cout << " == RESULTS FOR THE " << nIpiIterations << "-th iteration ==" << std::endl;
      std::cout.precision(10);
      maquis::cout << " - Energy difference for iteration = " << nIpiIterations << " = " << energyDifference << std::endl;
      maquis::cout << " - MPS overlap with solution at previous iteration = " << std::fabs(mpsOverlap) << std::endl;
      maquis::cout << std::endl;
      std::cout.precision(precision);
      // Checks convergence and, if not reached, starts a new IPI iteration
      if (nIpiIterations == numberOfOuterIterations || energyDifference < energyConvergenceThreshold ||
          std::fabs(1.-std::fabs(mpsOverlap)) < overlapConvergenceThreshold)
      {
        maquis::cout << " --> CONVERGENCE REACHED" << std::endl;
        convergedOuter = true;
      }
      else {
        maquis::cout << " --> CONVERGENCE NOT REACHED, STARTS NEW ITERATION" << std::endl;
        previousEnergy = nextEnergy;
        mpsBackup = this->mps;
      }
    }
  }

  /**
   * @brief Generic ALS-based optimization.
   * This routine solve a problem (that can be either a propagation, the solution of the linear system, or
   * the energy minimization) with the ALS algorithm. At the end, dumps the MPS, the energy, as well as the
   * measurements that were requested to be done at each microiteration.
   */
  void runAlternatingLeastSquares(std::string simulationType, int nSweeps, double energyThreshold)
  {
    // Reads in input parameters
    int meas_each = parms["measure_each"];
    int chkp_each = parms["chkp_each"];
    // -- Optimizer initialization --
    //std::shared_ptr<opt_base_t> optimizer;
    if (parms["optimization"] == "singlesite")
      // optimizer.reset( new ss_optimize<Matrix, SymmGroup, storage::disk>
      //                 (mps, mpo, parms, stop_callback, lat, init_site) );
      factory_ = std::make_unique<FactoryType>(simulationType, SweepOptimizationType::SingleSite, mps, mpo, parms, model, base::lat);
    else if(parms["optimization"] == "twosite")
      // optimizer.reset( new ts_optimize<Matrix, SymmGroup, storage::disk>
      //                 (mps, mpo, parms, stop_callback, lat, init_site) );
      factory_ = std::make_unique<FactoryType>(simulationType, SweepOptimizationType::TwoSite, mps, mpo, parms, model, base::lat);
    else
        throw std::runtime_error("Don't know this optimizer");
    // Retrieve the measurements that should be always done.
    auto always_measurements = this->iteration_measurements(init_sweep);
    auto firstEnergy = this->get_energy();
    energies_.push_back(firstEnergy);
    // Run the sweep-based simulation.
    try {
      for (int sweep=init_sweep; sweep < nSweeps; ++sweep) {
        // optimizer->sweep(sweep, Both);
        factory_->runSingleSweep(sweep);
        storage::disk::sync();
        bool converged = false;
        if ((sweep+1) % meas_each == 0 || (sweep+1) == nSweeps) {
          dumpParametersAndIterResults(sweep);
          dumpEnergy(sweep);
          if (!rfile().empty() && always_measurements.size() > 0)
            this->measure(this->results_archive_path(sweep) + "/results/", always_measurements);
          // stop simulation if an energy threshold has been specified
          int prev_sweep = sweep - meas_each;
          if (prev_sweep >= 0)
            converged = checkEnergyConvergence(energyThreshold);
        }
        last_sweep_ = sweep;
        /// write checkpoint
        bool stopped = stop_callback() || converged;
        if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == nSweeps)
          checkpoint_simulation(mps, sweep, -1);
        if (stopped)
          break;
      }
    }
    catch (dmrg::time_limit const& e) {
      maquis::cout << e.what() << " checkpointing partial result." << std::endl;
      checkpoint_simulation(mps, e.sweep(), e.site());
      dumpParametersAndIterResults(e.sweep());
      dumpEnergy(e.sweep());
    }
  }

  /** @brief Runs a propagation calculation */
  //AB For now it's mostly copy-pasted from optimize, should be rewritten in a cleaner way.
  /*
  void evolve()
  {
#ifdef DMRG_TD
    // Types definition
    using EvolverType = TimeEvolutionSweep<Matrix, SymmGroup, storage::disk>;
    using SSEvolverType = SingleSiteTimeEvolution<Matrix, SymmGroup, storage::disk>;
    using TSEvolverType = TwoSiteTimeEvolution<Matrix, SymmGroup, storage::disk>;
    // Reads in input parameters
    int meas_each = parms["measure_each"];
    int chkp_each = parms["chkp_each"];
    // Optimizer initialization
    std::shared_ptr<EvolverType> evolver;
    if (parms["optimization"] == "singlesite")
        evolver = std::make_shared<SSEvolverType>(mps, mpo, parms, stop_callback, init_site);
    else if(parms["optimization"] == "twosite")
        evolver = std::make_shared<TSEvolverType>(mps, mpo, parms, stop_callback, init_site);
    else
        throw std::runtime_error("Evolution modality not recognized");
    measurements_type always_measurements = this->iteration_measurements(init_sweep);
    int nSweeps = parms["nsweeps"];
    try {
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
            measurements_type always_measure = this->iteration_measurements(sweep);
            if (!parms["ALWAYS_MEASURE"].empty())
                this->measure(this->results_archive_path(sweep) + "/results/", always_measure);
            // Write iteration results
            {
              storage::archive ar(rfile(), "w");
              ar[this->results_archive_path(sweep) + "/results"] << evolver->iteration_results();
              ar[this->results_archive_path(sweep) + "/results/Energy/mean/value"] << std::vector<double>(1, energy);

            }
          }
        }
        // Dump data in
        last_sweep_ = sweep;
        bool stopped = stop_callback();
        if (stopped || (sweep + 1) % chkp_each == 0 || (sweep + 1) == parms["nsweeps"])
          checkpoint_simulation(mps, sweep, -1);
        if (stopped)
            break;
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
  void run_measure()
  {
    //if (this->get_last_sweep() < 0)
    //    throw std::runtime_error("Tried to measure before a sweep");
    this->measure("/spectrum/results/", all_measurements);
    // MPO creation
    MPO<Matrix, SymmGroup> mpoc = mpo;
    if (parms["use_compressed"])
        mpoc.compress(1e-12);
    double energy;
    // Measures the energy
    if (parms["MEASURE[Energy]"])
    {
        energy = maquis::real(expval(mps, mpoc))/maquis::real(overlap(mps, mps));
        maquis::cout << "Energy: " << energy << std::endl;
        maquis::cout << "MPS norm: " << maquis::real(overlap(mps, mps)) << std::endl;
        if (!rfile().empty()) {
            storage::archive ar(rfile(), "w");
            ar["/spectrum/results/Energy/mean/value"] << std::vector<double>(1, energy);
        }
    }
    // Measures the energy variance
    if (parms["MEASURE[EnergyVariance]"] > 0)
    {
        if (!parms["MEASURE[Energy]"])
            energy = maquis::real(expval(mps, mpoc));
        auto traitClass = MPOTimesMPSTraitClass<Matrix, SymmGroup>(mps, model, base::lat, model.total_quantum_numbers(parms),
                                                                   parms["max_bond_dimension"]);
        auto outputMPS = traitClass.applyMPO(mpoc);
        auto energy2 = maquis::real(overlap(outputMPS, outputMPS)/norm(mps));
        maquis::cout << "Energy^2: " << energy2 << std::endl;
        maquis::cout << "Variance: " << energy2 - energy*energy << std::endl;
        if (!rfile().empty())
        {
            storage::archive ar(rfile(), "w");
            ar["/spectrum/results/Energy^2/mean/value"] << std::vector<double>(1, energy2);
            ar["/spectrum/results/EnergyVariance/mean/value"] << std::vector<double>(1, energy2 - energy*energy);
        }
    }
    #if defined(HAVE_TwoU1) || defined(HAVE_TwoU1PG)
    if (!rfile().empty()) {
        BaseParameters parms_meas;
        parms_meas = parms.twou1_measurements();
        if (!parms_meas.empty())
            measure_transform<Matrix, SymmGroup>()(rfile(), "/spectrum/results", base::lat, mps, parms_meas);
    }
    else
        throw std::runtime_error("Transformed measurements not implemented yet without checkpoints");
    #endif
  }

  results_map_type measure_out() {
    results_map_type ret;
    // Do not measure before a sweep
    if (this->get_last_sweep() < 0)
      throw std::runtime_error("Tried to measure before a sweep");
    // Run all measurements and fill the result map
    for (auto&& meas: all_measurements)
      ret[meas.name()] = measure_and_save<Matrix,SymmGroup>(rfile(), "/spectrum/results", mps).meas_out(meas);
    // Measurements that require SU2U1->2U1 transformation
#if defined(HAVE_TwoU1) || defined(HAVE_TwoU1PG)
    BaseParameters parms_meas;
    parms_meas = parms.twou1_measurements();
    if (!parms_meas.empty()) {
      // Obtain a map with transformed measurements
      results_map_type transformed_meas = measure_transform<Matrix, SymmGroup>().meas_out(base::lat, mps, parms_meas, rfile(), "/spectrum/results");
      // Merge transformed measurements with the remaining results
      ret.insert(transformed_meas.begin(), transformed_meas.end());
    }
#endif
    return ret;
  }

  /** @brief Gets the energy for the mps that is stored in the sim object */
  /** if it is a feastMPS, return the feast energy of the zeroth feast state*/
  RealType get_energy() {
    if (!feastMPSs_)
      return maquis::real(expval(mps, mpo)/overlap(mps, mps));
    else
      return getFEASTEnergy(0);
  }

  /** @brief Gets the FEAST eigenstates - throws an exception if FEAST is not run */
  auto getFEASTEigenstates() {
    if (!feastMPSs_)
      throw std::runtime_error("FEAST eigenstate requested before running a FEAST simulation");
    else
      return feastMPSs_;
  }

  /** @brief Gets the FEAST energies -- throws an exception if FEAST is not run */
  RealType getFEASTEnergy(int iState) const {
    if (!feastMPSs_)
      throw std::runtime_error("FEAST energy requested before running a FEAST simulation");
    else if (iState >= feastMPSs_->size()) {
      std::string errorMessage = "FEAST energy requested for the"+std::to_string(iState)+"-th state, but only "+std::to_string(feastMPSs_->size())+" states are available";
      throw std::runtime_error(errorMessage);
    }
    else
      return maquis::real(expval(feastMPSs_->operator[](iState), mpo)/norm(feastMPSs_->operator[](iState)));
  }

  /**
   * @brief Method to extract a CI coefficient associated to a given determinant.
   *
   * Note that the format in which the input must be given is the same as for
   * the [hf_occ] (for electronic problems) or the [basis_state_generic] initializer
   * for the vibrational case.
   *
   * @param determinantString string associated with the target determinant
   * @return overlap value
   */
  typename Matrix::value_type getCICoefficient(std::string determinantString) override {
      auto modifiedParameters = parms;
      std::string initState = (parms["MODEL"] == "quantum_chemistry") ? "hf" : "basis_state_generic";
      modifiedParameters.set("init_type", initState);
      if (parms["MODEL"] == "quantum_chemistry")
          modifiedParameters.set("hf_occ", determinantString);
      else
          modifiedParameters.set("init_basis_state", determinantString);
      auto mpsOverlap = MPSType(lat.size(), *(model.initializer(lat, modifiedParameters)));
      return overlap(mpsOverlap, mps)/std::sqrt(norm(mpsOverlap)*norm(mps));
  }

  /** @brief Updates the integral and regenerates the data that depends on it */
  void update_integrals(const chem::integral_map<typename Matrix::value_type> & integrals)
  {
      if (parms.is_set("integral_file") || parms.is_set("integrals"))
          throw std::runtime_error("updating integrals in the interface not supported yet in the FCIDUMP format");
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

  results_collector& get_iteration_results()
  {
    // If iteration_results is empty, we didn't perform the sweep yet, but possibly loaded the MPS from a checkpoint
    // so we need to load also iteration results
    if (iteration_results_.empty())
    {
      // If we are not loading from a checkpoint, last_sweep is set to -1
      // so we need to return an empty iteration_results vector
      if (get_last_sweep() < 0)
        return iteration_results_;
      // otherwise, we are restarting but there's something wrong with the checkpoint
      if (!rfile().empty()) {
        try { // Load the iteration results from the last sweep
          storage::archive ar(rfile(), "r");
          ar[results_archive_path(last_sweep_) + "/results"] >> iteration_results_;
        }
        catch (std::exception& e) {
          maquis::cerr << e.what() << std::endl;
          throw std::runtime_error("Error reading iteration results from checkpoint.");
        }
      }
      else
          throw std::runtime_error("No result file specified for restart -- cannot read iteration results!");
    }
    return iteration_results_;
  }

  /** @brief Get the overlap of the MPS with another MPS, which is loaded from a chkp file */
  virtual typename Matrix::value_type get_overlap(const std::string & aux_filename)
  {
      maquis::checks::symmetry_check(parms, aux_filename);
      MPS<Matrix, SymmGroup> aux_mps;
      load(aux_filename, aux_mps);
      return overlap(aux_mps, this->mps);
  }

  /** @brief Getter for the number of sweeps that have been run */
  int get_last_sweep() { return last_sweep_; };

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
    if (!rfile().empty()) {
      auto energy = this->get_energy();
      energies_.push_back(energy);
      storage::archive ar(rfile(), "w");
      ar[this->results_archive_path(iSweep) + "/results/Energy/mean/value"] << std::vector<double>(1, energy);
    }
  }

  /**  @brief Checks energy convergence of the sweep-based optimization */
  bool checkEnergyConvergence(double convergenceThreshold) {
    bool converged = false;
    auto emin = *std::min_element(energies_.begin(), energies_.end()-1);
    auto eminNew = *std::min_element(energies_.begin(), energies_.end());
    auto eDiff = std::abs(emin - eminNew);
    if (eDiff < convergenceThreshold)
      converged = true;
    return converged;
  }

  /** @brief Returns the path where the result of a given sweep are stored */
  std::string results_archive_path(int sweep) const {
    status_type status;
    status["sweep"] = sweep;
    return base::results_archive_path(status);
  }

  /** @brief Dumps the simulation results to the checkpoint file */
  void checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, int sweep, int site, std::string filename = "") {
    status_type status;
    status["sweep"] = sweep;
    status["site"]  = site;
    return base::checkpoint_simulation(state, status, filename);
  }

  void dumpParameters(std::string filename = "") {
    if (!chkpfolder().empty()) {
      std::string chkpfilename;
      if (filename.empty())
        chkpfilename = chkpfolder();
      else
        chkpfilename = chkpfolder() + "_" + filename;
      storage::archive ar(chkpfilename+"/props.h5", "w");
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
