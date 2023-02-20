/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SWEEP_BASED_LINEAR_SYSTEM_H
#define SWEEP_BASED_LINEAR_SYSTEM_H

#include "GenericSweepSimulation.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/LinearSystem/linsolver.h"
#include "dmrg/LinearSystem/LinSystemTraitsClass.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "dmrg/utils/checks.h"
#include "BoundaryPropagator.h"
#include "OverlapPropagator.h"
#include "SweepOptimizationTypeTrait.h"

template<class Matrix, class SymmGroup, class Storage, SweepOptimizationType SweepType>
class SweepBasedLinearSystem : public GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType> {
public:
  using Base = GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType>;
  using OverlapPropagatorType = OverlapPropagator<Matrix, SymmGroup, Storage>;
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepType>;
  using SiteProblemType = SiteProblem<Matrix, SymmGroup>;
  using LinearSolverType = LinSolver<Matrix, SymmGroup>;
  using ModelType =  typename Base::ModelType;
  using MPSType = typename Base::MPSType;
  using MPOType = typename Base::MPOType;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using BlockMatrixType = block_matrix<Matrix, SymmGroup>;
  using ValueType = typename MPSTensorType::scalar_type;
  //
  using Base::boundaryPropagator_;
  using Base::indexOfMicroIteration_;
  using Base::iterationResults_;
  using Base::lattice_;
  using Base::L_;
  using Base::model_;
  using Base::mps_;
  using Base::mpsContainer_;
  using Base::mpoContainer_;
  using Base::parms_;
  using Base::siteLeft_;
  using Base::siteRight_;
  using Base::verbose_;

  /** @brief Class constructor */
  SweepBasedLinearSystem(MPSType& mps, const MPOType& mpo, BaseParameters& parms, const ModelType& model,
                         const Lattice& lattice, bool verbose)
    : Base(mps, mpo, parms, model, lattice, verbose, std::string("Linear system solver")),
      adaptiveBondDimension_(false), shiftParameter_(0.), isPrecond_(false), rhsMps_(mps)
  {
    /* // Folded simulation --> To be reactivated when implementing the folded operator
    if (parms["pI_folded"] == "yes") {
        maquis::cout << " Activating folded treatment " << std::endl;
        isSquared = true;
    } */
    overlapPropagator_ = std::make_unique<OverlapPropagatorType>(mps_, rhsMps_);
    /* To be reactivated when implementing the folded operator
    if (isSquared) {
      leftSquared_.resize(mpo.length()+1);
      rightSquared_.resize(mpo.length()+1);
      leftCross_.resize(mpo.length()+1);
      rightCross_.resize(mpo.length()+1);
    } */
    // Adaptive m
    if (parms.is_set("linsystem_truncation_ratio")) {
      adaptiveBondDimension_ = true;
      truncationRatio_ = parms["linsystem_truncation_ratio"].as<double>();
    }
    // Note that we subtract the core energy to the shift parameter (the SiteProblem object
    // does not include that contribution)
    if (parms_.is_set("ipi_shift"))
      shiftParameter_ = parms["ipi_shift"].as<ValueType>()-mpoContainer_.getMPO().getCoreEnergy();
    if (parms_["linsystem_precond"] == "yes")
      isPrecond_ = true;
    calculateExactError_ = (parms_["linsystem_exact_error"] == "yes");
  }

  /** @brief Setter for the shift */
  void setShift(ValueType newShift) {
    shiftParameter_ = newShift;
  }

  /** @brief Method called at the beginning of each sweep */
  void prepareSweep() override final {
    iterationResults_.clear();
  }

  /** @brief Method called before each microiteration */
  void prepareMicroiteration() override final {
    siteProblem_ = std::make_unique<SiteProblemType>(boundaryPropagator_->getLeftBoundary(siteLeft_), boundaryPropagator_->getRightBoundary(siteRight_),
                                                     mpoContainer_.getMPOTensor(siteLeft_));
    rhs_ = overlapPropagator_->template getOrthogonalVector<SweepType>(siteLeft_, siteRight_);
    if (isPrecond_)
      preconditioner_ = std::make_unique<BlockMatrixType>(contraction::diagonal_hamiltonian(boundaryPropagator_->getLeftBoundary(siteLeft_),
                                                                                            boundaryPropagator_->getRightBoundary(siteRight_),
                                                                                            mpoContainer_.getMPOTensor(siteLeft_),
                                                                                            mpsContainer_.getMPSTensor(siteLeft_)));
  }

  /** @brief Solution of the site-centered problem */
  MPSTensorType solveLocalProblem() override final {
    auto& mpsToOptimize = mpsContainer_.getMPSTensor(siteLeft_);
    LinearSolverType ls(siteProblem_, mpsToOptimize, rhs_, shiftParameter_, parms_, preconditioner_, verbose_);
    auto resultOfLocalSiteProblem = ls.res();
    iterationResults_["Energy"] << std::get<0>(resultOfLocalSiteProblem) + maquis::real(mpoContainer_.getMPO().getCoreEnergy());
    energyPerMicroIter_.push_back(std::get<0>(resultOfLocalSiteProblem));
    errorPerMicroIter_.push_back(std::get<1>(resultOfLocalSiteProblem));
    return std::get<2>(resultOfLocalSiteProblem);
  }

  /** @brief Propagates the orthogonal vector */
  void propagateOtherTensors() override final {
    auto sweepType = SweepTraitClass::getSweepDirection(L_, indexOfMicroIteration_);
    // Boundary propagation
    if (sweepType == SweepDirectionType::Forward &&
        !SweepTraitClass::changeDirectionNextMicroiteration(L_, indexOfMicroIteration_)) {
      rhsMps_.move_normalization_l2r(siteLeft_, siteLeft_+1);
      if (overlapPropagator_)
        overlapPropagator_->updateLeftOverlapBoundaries(siteLeft_+1);
    }
    else {
      auto mpsCopy = rhsMps_;
      rhsMps_.move_normalization_r2l(siteRight_-1, siteRight_-2);
      if (overlapPropagator_)
        overlapPropagator_->updateRightOverlapBoundaries(siteRight_-1);
    }
  }

  /** @brief Operations to be executed at the end of a microiteration */
  void finalizeMicroIteration(const truncation_results& trunc) override final {
    iterationResults_["BondDimension"]   << trunc.bond_dimension;
    iterationResults_["TruncatedWeight"] << trunc.truncated_weight;
    iterationResults_["SmallestEV"]      << trunc.smallest_ev;
  }

  /** @brief Operations to be executed at the end of the sweep */
  void finalizeSweep() override final {
    if (calculateExactError_) {
      int mMax = parms_["max_bond_dimension"];
      auto error = LinSystemTraitClass<Matrix, SymmGroup>::calculateError(mpsContainer_.getMPS(), rhsMps_, mpoContainer_.getMPO(), shiftParameter_,
                                                                          model_, lattice_, model_.total_quantum_numbers(parms_), mMax);
      if (verbose_) {
        maquis::cout << std::scientific << std::setprecision(16);
        maquis::cout << " Exact error = " << error << std::endl;
        maquis::cout << std::endl;
      }
    }
  }

  /** @brief Prints a summary of the results of a minimization */
  void printSummary() const {
    // Prints header
    maquis::cout << std::endl;
    maquis::cout << " == SUMMARY OF THE SWEEP-BASED SOLUTION OF THE LINEAR SYSTEM == " << std::endl;
    maquis::cout << std::endl;
    maquis::cout << " +----------------+------------------+------------------+" << std::endl;
    maquis::cout << "   Microiteration |      Energy      |       Error       " << std::endl;
    maquis::cout << " +----------------+------------------+------------------+" << std::endl;
    for (int iIter = 0; iIter < energyPerMicroIter_.size(); iIter++) {
      maquis::cout << std::setw(17) << std::right << iIter
                   << std::setw(19) << std::setprecision(10) << std::scientific << std::right << energyPerMicroIter_[iIter]
                   << std::setw(19) << std::setprecision(10) << std::scientific << std::right << errorPerMicroIter_[iIter] << std::endl;
      if ((iIter+1)%(SweepTraitClass::getNumberOfMicroiterations(L_)) == 0)
        maquis::cout << " +----------------+------------------+------------------+" << std::endl;
    }
    maquis::cout << std::endl;
  }

  /** @brief Whether to normalize the MPS at the end of a half-sweep */
  bool normalizeAtEnd() override final {
    return false;
  }

  /** @brief Gets the system rhs */
  auto getRhs() const {
    return rhsMps_;
  }

private:
  // Class members
  MPSType rhsMps_;                                                // RHS for the solution of the linear system.
  std::shared_ptr<BlockMatrixType> preconditioner_;               // If needed, stores the preconditioner.
  bool adaptiveBondDimension_;                                    // Whether to dynamically adapt the bond dimension.
  bool isPrecond_;                                                // If true, activates the preconditioning.
  bool calculateExactError_;                                      // If true, calculates the exact error associated to the solution of the linear system.
  double truncationRatio_;                                        // Parameter for a DBSS-like solution of the linear system.
  ValueType shiftParameter_;                                      // Shift parameter for the linear system
  MPSTensorType rhs_;                                             // RHS of the local linear system (updated at each microiteration).
  std::unique_ptr<OverlapPropagatorType> overlapPropagator_;      // Object needed to store the partial MPS/MPS contraction
  std::shared_ptr<SiteProblemType> siteProblem_;                  // Site problem associated with the solution of the linear system.
  std::vector<double> energyPerMicroIter_, errorPerMicroIter_;    // Backup of results along the propagation.
};

#endif // SWEEP_BASED_LINEAR_SYSTEM_H
