/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SWEEP_BASED_ENERGY_MINIMIZATION_H
#define SWEEP_BASED_ENERGY_MINIMIZATION_H

#include "GenericSweepSimulation.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/optimize/ietl_jacobi_davidson.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "dmrg/utils/checks.h"
#include "BoundaryPropagator.h"
#include "OverlapPropagator.h"

template<class Matrix, class SymmGroup, class Storage, SweepOptimizationType SweepType>
class SweepBasedEnergyMinimization : public GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType> {
public:
  using Base = GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType>;
  using OverlapPropagatorType = OverlapPropagator<Matrix, SymmGroup, Storage>;
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepType>;
  using SiteProblemType = SiteProblem<Matrix, SymmGroup>;
  using MPSType = typename Base::MPSType;
  using ModelType = typename Base::ModelType;
  using MPOType = typename Base::MPOType;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using ValueType = typename MPSTensorType::value_type;
  //
  using Base::boundaryPropagator_;
  using Base::getSpecificResult;
  using Base::indexOfMicroIteration_;
  using Base::iterationResults_;
  using Base::L_;
  using Base::mps_;
  using Base::mpsContainer_;
  using Base::mpoContainer_;
  using Base::parms_;
  using Base::runSweepSimulation;
  using Base::siteLeft_;
  using Base::siteRight_;

  /** @brief Class constructor */
  SweepBasedEnergyMinimization(MPSType& mps, const MPOType& mpo, BaseParameters& parms, const ModelType& model,
                               const Lattice& lattice, bool verbose)
    : Base(mps, mpo, parms, model, lattice, verbose, std::string("Optimization")), nOrtho_(0)
  {
    if (parms_.is_set("ortho_states") && parms_["ortho_states"] != "") {
      files_ = parms_["ortho_states"].str();
      std::vector<std::string> files;
      boost::split(files, files_, boost::is_any_of(", "));
      if (!parms_.is_set("n_ortho_states"))
        throw std::runtime_error("Please set [n_ortho_states]");
      else
        nOrtho_ = parms_["n_ortho_states"];
      overlapPropagator_ = std::make_unique<OverlapPropagatorType>(mps_, files, parms_);
      if (nOrtho_ != overlapPropagator_->getNumberOfOverlapMPSs())
        throw std::runtime_error("Nuber of chkp files not coherent with [n_ortho_states] parameter");
      orthoLocal_.resize(nOrtho_);
      maquis::cout << "Running a constrained optimization with respect to " << nOrtho_ << " states." << std::endl;
    }
  }

  /** @brief Method called at the beginning of each sweep */
  void prepareSweep() override final {
    iterationResults_.clear();
  }

  /** @brief Method called before each microiteration */
  void prepareMicroiteration() override final {
    siteProblem_ = std::make_unique<SiteProblemType>(boundaryPropagator_->getLeftBoundary(siteLeft_), boundaryPropagator_->getRightBoundary(siteRight_),
                                                     mpoContainer_.getMPOTensor(siteLeft_));
    // std::cout << "Initial energy " << siteProblem_->get_energy(mpsContainer_.getMPSTensor(siteLeft_)) + mpoContainer_.getMPO().getCoreEnergy() << std::endl;
    if (overlapPropagator_)
      for (int iState = 0; iState < nOrtho_; iState++)
        orthoLocal_[iState] = overlapPropagator_->template getOrthogonalVector<SweepType>(iState, siteLeft_, siteRight_);
  }

  /** @brief Solution of the site-centered problem */
  MPSTensorType solveLocalProblem() override final {
    auto& mpsToOptimize = mpsContainer_.getMPSTensor(siteLeft_);
    if (parms_["eigensolver"] == std::string("IETL"))
      resultOfLocalSiteProblem_ = solve_ietl_lanczos(*(siteProblem_.get()), mpsToOptimize, parms_);
    else if (parms_["eigensolver"] == std::string("IETL_JCD"))
      resultOfLocalSiteProblem_ = solve_ietl_jcd(*(siteProblem_.get()), mpsToOptimize, parms_, orthoLocal_);
    else if (parms_["eigensolver"] == std::string("IETL_DAVIDSON"))
      resultOfLocalSiteProblem_ = solve_ietl_jcd(*(siteProblem_.get()), mpsToOptimize, parms_, orthoLocal_);
    else
      throw std::runtime_error("I don't know this eigensolver.");
    // Loads the final results
    auto energy = resultOfLocalSiteProblem_.first + mpoContainer_.getMPO().getCoreEnergy();
    maquis::cout << std::setprecision(10) << " Energy = " << std::setprecision(16) << energy << std::endl;
    iterationResults_["Energy"] << energy;
    return resultOfLocalSiteProblem_.second;
  }

  /** @brief Propagates the boundaries */
  void propagateOtherTensors() override final {
    auto sweepType = SweepTraitClass::getSweepDirection(L_, indexOfMicroIteration_);
    // Boundary propagation
    if (overlapPropagator_) {
      if (sweepType == SweepDirectionType::Forward &&
          !SweepTraitClass::changeDirectionNextMicroiteration(L_, indexOfMicroIteration_)) {
        overlapPropagator_->updateLeftOverlapBoundaries(siteLeft_+1);
      }
      else {
        overlapPropagator_->updateRightOverlapBoundaries(siteRight_-1);
      }
    }
  }

  /** @brief Operations to be executed at the end of a microiteration */
  void finalizeMicroIteration(const truncation_results& trunc) override final {
    iterationResults_["BondDimension"]   << trunc.bond_dimension;
    iterationResults_["TruncatedWeight"] << trunc.truncated_weight;
    iterationResults_["SmallestEV"]      << trunc.smallest_ev;
  }

  /** @brief Operations to be executed at the end of the sweep */
  void finalizeSweep() override final { }

  /** @brief Whether to normalize the MPS at the end of a half-sweep */
  bool normalizeAtEnd() override final {
    return true;
  }

private:
  // Class members
  int nOrtho_;
  std::vector<MPSTensorType> orthoLocal_;
  std::unique_ptr<OverlapPropagatorType> overlapPropagator_;
  std::unique_ptr<SiteProblemType> siteProblem_;
  std::string files_;
  std::pair<ValueType, MPSTensorType > resultOfLocalSiteProblem_;
};

#endif // SWEEP_BASED_ENERGY_MINIMIZATION_H
