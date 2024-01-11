/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef GENERIC_SWEEPS_SIMULATION_H
#define GENERIC_SWEEPS_SIMULATION_H

#include <any>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "dmrg/SweepBasedAlgorithms/BoundaryPropagator.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/results_collector.h"
#include "SweepMPSContainer.h"
#include "SweepMPOContainer.h"
#include "SweepMPSUpdater.h"
#include "SweepOptimizationTypeTrait.h"
#include "utils/io.hpp"

/**
 * @brief Class representing a generic sweep-based simulation.
 *
 * @tparam Matrix matrix class underlying the matrix storage.
 * @tparam SymmGroup symmetry group of the Hamiltonian.
 */

template<class Matrix, class SymmGroup, class Storage, SweepOptimizationType SweepType>
class GenericSweepSimulation {
public:
  // Types declaration
  using ModelType = Model<Matrix, SymmGroup>;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using MPSContainerType = SweepMPSContainer<Matrix, SymmGroup, SweepType>;
  using MPOContainerType = SweepMPOContainer<Matrix, SymmGroup, SweepType>;
  using SweepMPSUpdaterType = SweepMPSUpdater<Matrix, SymmGroup, Storage, SweepType>;
  using BoundaryPropagatorType = BoundaryPropagator<Matrix, SymmGroup, Storage>;
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepType>;

  /** @brief Class constructor */
  GenericSweepSimulation(MPSType& mps, const MPOType& mpo, BaseParameters& parms, const ModelType& model,
                         const Lattice& lattice, bool verbose, std::string simulationName="Optimization")
    : mps_(mps), mpoContainer_(mpo, mps), mpsContainer_(mps),
      L_(mps_.length()), indexOfMicroIteration_(0), siteLeft_(0), siteRight_(1),
      nSweeps_(0), parms_(parms), simulationName_(std::move(simulationName)),
      model_(model), lattice_(lattice), verbose_(verbose)
  {
    
    
    mps_.normalize_right();
    nSweeps_ = parms_["nsweeps"];
    boundaryPropagator_ = std::make_shared<BoundaryPropagatorType>(mps_, mpoContainer_.getMPO());
    
    mpsUpdater_ = std::make_unique<SweepMPSUpdaterType>(mpoContainer_.getMPO(), mps_, boundaryPropagator_, parms_, verbose_);
  }

  /** @brief Virtual destructor */
  virtual ~GenericSweepSimulation() = default; 
  /**
   * @brief Execution of a generic sweep-based optimization algorithm.
   *
   * Note that we delegate every action to the derived class, with the exception of the
   * memory management, which is done here to ensure that the implementation is consistent
   * for all methods.
   */
  void runSweepSimulation() {
    // == LOOP OVER THE SWEEPS ==
    for (int iSweep = 0; iSweep < nSweeps_; iSweep++) {
      this->runSingleSweep(iSweep);
    }
  }


  /** @brief Runs a single sweep of a sweep-based optimization */
  void runSingleSweep(int iSweep) {
    // Prints header
    if (iSweep == 0) { printGenericInfo(); }
    // Preparatory operations.
    this->prepareSweep();
    indexOfMicroIteration_ = 0;
    this->printSweepSpecificInfo(iSweep);
    this->updateSites(); // Needed for initialization
    // Prefetches the boundaries that will be needed for the first sweep
    Storage::prefetch(boundaryPropagator_->getLeftBoundary(siteLeft_));
    Storage::prefetch(boundaryPropagator_->getRightBoundary(siteRight_));
    auto numberOfMicroIterations = SweepTraitClass::getNumberOfMicroiterations(L_);
    auto start_sweep = std::chrono::high_resolution_clock::now();
    // == LOOP OVER THE MICROITERATIONS ==
    while (indexOfMicroIteration_ < numberOfMicroIterations) {
      auto start = std::chrono::high_resolution_clock::now();
      // Useful local variables
      auto sweepType = SweepTraitClass::getSweepDirection(L_, indexOfMicroIteration_);
      auto changeDirection = SweepTraitClass::changeDirectionNextMicroiteration(L_, indexOfMicroIteration_);
      // Syncs the storage
      // if (sweepType == SweepDirectionType::Backward && indexOfMicroIteration_ == lastSite_) {
      //   std::cout << "Synchronizing storage" << std::endl;
      //   Storage::sync();
      // }
      this->updateSites();
      if (verbose_) {
        printMicroiterInfo(sweepType);
      }
      // Gets the boundary that are needed. Note that, in a forward sweep, the left boundary is assumed
      // to have been generated during the previous boundary update and, therefore, is not fetched.
      if (sweepType == SweepDirectionType::Backward || indexOfMicroIteration_ == 0) {
        Storage::fetch(boundaryPropagator_->getLeftBoundary(siteLeft_));
      }
      if (sweepType == SweepDirectionType::Forward) {
        Storage::fetch(boundaryPropagator_->getRightBoundary(siteRight_));
      }
      // Starts prefetching what will be needed in the following microiteration.
      // Note that, for instance, we don't prefetch the left boundary for the l2r sweep because
      // this will be taken care in the boundary propagation (in other words, there is no
      // need to prefetch the left boundary since it will be anyways modified by the boundary
      // propagation)
      if (sweepType == SweepDirectionType::Forward) {
        auto sweepTypeNext = SweepTraitClass::getSweepDirection(L_, indexOfMicroIteration_+1);
        if (sweepTypeNext == sweepType) {
          Storage::prefetch(boundaryPropagator_->getRightBoundary(SweepTraitClass::getIndexOfRightBoundary(L_, indexOfMicroIteration_+1)));
        } else if (sweepTypeNext != SweepDirectionType::EndOfLattice) {
          Storage::prefetch(boundaryPropagator_->getLeftBoundary(SweepTraitClass::getIndexOfLeftBoundary(L_, indexOfMicroIteration_+1)));
        }
      }
      else if (sweepType == SweepDirectionType::Backward) {
        if (indexOfMicroIteration_ != numberOfMicroIterations-1) {
          Storage::prefetch(boundaryPropagator_->getLeftBoundary(SweepTraitClass::getIndexOfLeftBoundary(L_, indexOfMicroIteration_+1)));
        }
      }
      // == SOLUTION OF THE LOCAL PROBLEM ==
      this->prepareMicroiteration();
      auto outputTensor = this->solveLocalProblem();
      // == MPS UPDATE ==
      auto boundaryGrowthModality = (sweepType == SweepDirectionType::Forward && !changeDirection) ? GrowBoundaryModality::LeftToRight
                                                                                                   : GrowBoundaryModality::RightToLeft;
      auto truncationResults = mpsUpdater_->generateUnitaryFactor(
          siteLeft_, siteRight_,
          boundaryGrowthModality,
          outputTensor,
          this->getAlpha(iSweep), this->get_cutoff(iSweep),
          this->get_Mmax(iSweep), this->normalizeAtEnd(),
          this->activatePerturbation());
      // == BOUNDARY PROPAGATION ==
      // Updates the boundary
      this->propagateBoundaries();
      this->propagateOtherTensors();
      this->performBackPropagation(boundaryGrowthModality);
      mpsUpdater_->mergeUnitaryFactor(boundaryGrowthModality, siteLeft_, siteRight_, this->normalizeAtEnd());
      // Now, can drop the memory of the right boundary (in the case of a l2r sweep).
      // The memory will anyways be overwritten by the r2l sweep that will follow.
      // Note also that, if we are at a point at which we reverse the direction of the boundary
      // propagation, we don't drop the right boundary because the next step will be a r2l sweep
      // and, therefore, that element of the right boundary won't be overwritten.
      if (sweepType == SweepDirectionType::Forward && !changeDirection) {
        Storage::drop(boundaryPropagator_->getRightBoundary(siteRight_));
      } else { // if (sweepType == SweepDirectionType::Backward)
        Storage::drop(boundaryPropagator_->getLeftBoundary(siteLeft_));
      }
      this->finalizeMicroIteration(truncationResults);
      indexOfMicroIteration_ += 1;
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
      if (verbose_) {
        maquis::cout << "[Time: " << std::setprecision(2) << duration.count()
          << std::setprecision(6) << " s]\n";
      }
    }
    // At the end, just stores to file the final right boundary (if needed, one can use it for the next sweep)
    // Storage::StoreToFile(boundaryPropagator_->getRightBoundary(siteRight_-1));

    auto stop_sweep = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duriation_sweep = stop_sweep - start_sweep;
    this->finalizeSweep();
    printSweepEnergy();
    maquis::cout << "[Sweep " << iSweep << " took "
      << std::setprecision(2) << duriation_sweep.count() << std::setprecision(6)
      << " s]\n\n";
  }

  /** @brief Gets the container with the results of each iteration */
  auto iteration_results() const { return iterationResults_; }

  /** @brief Gets a specific value of the iteration result */
  template<class CastType>
  CastType getSpecificResult(const std::string& resultName) {
    if (!iterationResults_.has(resultName)) {
      throw std::runtime_error("Trying to access non-existing simulation result");
    }
    return std::any_cast<CastType>(iterationResults_[resultName].get()[0]);
  }

  /** @brief Prints final energy and difference in energy between final
   *         and first microiteration */
  void printSweepEnergy() {
    if (!iterationResults_.has("Energy")) {
      throw std::runtime_error("Trying to access non-existing 'Energy' result");
    }
    std::vector<std::any> energy_vec = iterationResults_["Energy"].get();
    if (energy_vec.back().type() == typeid(double)) {
      auto energy = std::any_cast<double>(energy_vec.back());
      maquis::cout << " SWEEP ENERGY = " << std::setprecision(16) << energy << " ";
      if (energy_vec.size() > 1) {
        auto energy_prev = std::any_cast<double>(energy_vec.front());
        maquis::cout << std::scientific << std::setprecision(3);
        maquis::cout << "  (Delta E = " << energy - energy_prev << ")  ";
        maquis::cout << std::defaultfloat;
      }
    } else if (energy_vec.back().type() == typeid(std::complex<double>)) {
      auto energy = std::any_cast<std::complex<double>>(energy_vec.back());
      maquis::cout << " SWEEP ENERGY (Complex Valued) = " << std::setprecision(16) << std::real(energy) << " + i " << std::imag(energy) << " ";
      if (energy_vec.size() > 1) {
        auto energy_prev = std::any_cast<std::complex<double>>(energy_vec.front());
        double real_diff = std::real(energy) - std::real(energy_prev);
        double imag_diff = std::imag(energy) - std::imag(energy_prev);
        maquis::cout << std::scientific << std::setprecision(3);
        maquis::cout << "  (Delta E = " << real_diff << " + i " << imag_diff << ")  ";
        maquis::cout << std::defaultfloat;
      }
    }
  }

protected:

  /**
   * @brief Collects the operation to be done before a sweep.
   * Note that these operations are done only once per sweep, i.e. they are not
   * repeated at each microiteration.
   */
  virtual void prepareSweep() = 0;

  /** @brief Collects the operation to be performed before a microiteration */
  virtual void prepareMicroiteration() = 0;

  /** @brief Runs the actual sweep simulation */
  virtual MPSTensorType solveLocalProblem() = 0;

  /** @brief Whether the MPS should be normalized at the end of a sweep */
  virtual bool normalizeAtEnd() = 0;

  /**
   * @brief Whether to activate noise.
   *
   * By default, noise is activated. However, this choice can be funneled to the
   * derived class -- this is the case for TD simulations, where the noise may
   * compromise the energy conservation.
   */
  virtual bool activatePerturbation() { return true; }

  /** @brief Collects the operation to be run at the end of a micro iteration */
  virtual void finalizeMicroIteration(const truncation_results& trunc) = 0;

  /** @brief Collects the operation to be run at the end of the simulation */
  virtual void finalizeSweep() = 0;

  /**
   * @brief Propagates other tensor networks that may be needed.
   * Note that this method is called *after* the back-propagation.
   */
  virtual void propagateOtherTensors() {}

  /**
   * @brief Back-propagation step.
   * By default, nothing is done. However, this function can be defined by
   * derived class to perform a back-propagation step.
   */
  virtual void performBackPropagation(GrowBoundaryModality boundaryGrowthModality) {};

  /** @brief Boundary propagation method */
  void propagateBoundaries() {
    auto sweepType = SweepTraitClass::getSweepDirection(L_, indexOfMicroIteration_);
    // Boundary propagation
    if (sweepType == SweepDirectionType::Forward &&
        !SweepTraitClass::changeDirectionNextMicroiteration(L_, indexOfMicroIteration_)) {
      boundaryPropagator_->updateLeftBoundary(siteLeft_+1);
    } else {
      boundaryPropagator_->updateRightBoundary(siteRight_-1);
    }
  };

    /** @brief Simple utility function for a logarithmic interpolation */
  static double log_interpolate(double y0, double y1, int N, int i)
  {
    if (y0 <= 0) // Safeguard if for example -1 is entered
      return 0;
    if (i == 0)
      return y0;
    if (i >= N)
      return y1;
    if (y1 <= 0) // Safeguard if for example -1 is entered
      y1 = 1e-16;
    double x = log(y1/y0)/N;
    return y0*exp(x*i);
  }


  /** @brief Method to get the truncation threshold for a given sweep */
  double get_cutoff(int sweep) const {
    return (sweep >= parms_.template get<int>("ngrowsweeps")) ? parms_.template get<double>("truncation_main")
        : log_interpolate(parms_.template get<double>("truncation_initial"),
                          parms_.template get<double>("truncation_main"),
                          parms_.template get<int>("ngrowsweeps"), sweep);
  }

  /** @brief Method to get the proper noise parameter */
  auto getAlpha(int iSweep) const {
    double alpha;
    int ngs = parms_.template get<int>("ngrowsweeps");
    int nms = parms_.template get<int>("nmainsweeps");
    if (iSweep < ngs) {
      alpha = parms_.template get<double>("alpha_initial");
    } else if (iSweep < ngs + nms) {
      alpha = parms_.template get<double>("alpha_main");
    } else {
      alpha = parms_.template get<double>("alpha_final");
    }
    return alpha;
  }

  /** @brief Calculats the maximum value of the bond dimension for a given sweep */
  int get_Mmax(int sweep) const {
    std::size_t Mmax;
    if (parms_.is_set("sweep_bond_dimensions")) {
      std::vector<std::size_t> ssizes = parms_.template get<std::vector<std::size_t> >("sweep_bond_dimensions");
      if (sweep >= ssizes.size()) {
        Mmax = *ssizes.rbegin();
      } else {
        Mmax = ssizes[sweep];
      }
    } else {
      Mmax = parms_.template get<std::size_t>("max_bond_dimension");
    }
    return Mmax;
  }

  /** @brief Updates the index of the sites */
  void updateSites() {
    siteLeft_ = SweepTraitClass::getIndexOfLeftBoundary(L_, indexOfMicroIteration_);
    siteRight_ = SweepTraitClass::getIndexOfRightBoundary(L_, indexOfMicroIteration_);
  }

  /** @brief Prints generic information about the */
  void printGenericInfo() const {
    maquis::cout << "+--------------------------------------------------+\n"
                 << " " << SweepTraitClass::getSimulationTypeName() 
                 << " " << simulationName_ << '\n'
                 << " - Max Sweeps: " << nSweeps_ << '\n'
                 << "+--------------------------------------------------+\n";
  }

  /** @brief Prints info that are sweep-specific */
  void printSweepSpecificInfo(int iSweep) const {
    maquis::cout << std::defaultfloat;
    maquis::cout << "SWEEP " << iSweep + 1 << "; Noise: " << this->getAlpha(iSweep)
                 << "; Max Bond Dim: " << this->get_Mmax(iSweep) 
                 << "; Truncation: " << this->get_cutoff(iSweep) << '\n';
  }

  /** @brief Prints information regarding the current microiteration */
  void printMicroiterInfo(SweepDirectionType sweepType) const {
    maquis::cout << std::flush;
    maquis::cout << (sweepType == SweepDirectionType::Backward ? "  <" : "   ") << "-" << siteLeft_ << "-";
    if (SweepType == SweepOptimizationType::TwoSite) {
      maquis::cout << (siteLeft_ + 1) << "-"; 
    }
    maquis::cout << (sweepType == SweepDirectionType::Forward ? "> " : "  ");
  }

  /** @brief Checks whether the current microiteration is associated with a terminal site */
  inline bool isTerminal() const {
    return indexOfMicroIteration_ == 0 || (SweepTraitClass::changeDirectionNextMicroiteration(L_, indexOfMicroIteration_)
                                          && SweepTraitClass::getSweepDirection(L_, indexOfMicroIteration_) == SweepDirectionType::Forward);
  }

protected:
  MPSType& mps_;
  MPOContainerType mpoContainer_;
  MPSContainerType mpsContainer_;
  std::unique_ptr<SweepMPSUpdaterType> mpsUpdater_;
  int L_, indexOfMicroIteration_, siteLeft_, siteRight_, nSweeps_;
  BaseParameters& parms_;
  results_collector iterationResults_;
  std::shared_ptr<BoundaryPropagatorType> boundaryPropagator_;
  std::string simulationName_;
  const ModelType& model_;
  const Lattice& lattice_;
  bool verbose_;
};

#endif // GENERIC_SWEEPS_SIMULATION_H
