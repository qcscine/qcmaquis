/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef FEAST_SIMULATOR
#define FEAST_SIMULATOR

#include <cstdlib>
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/SweepBasedAlgorithms/SweepBasedLinearSystem.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/storage.h"
#include "FEASTQuadrature.h"
#include "FEASTPostProcessor.h"

template<class SymmGroup>
class FEASTSimulator {
  using Matrix = cmatrix;
  using StorageType = storage::disk; // This is hardcoded for now
  using ComplexType = std::complex<double>;
  using LatticeType = Lattice;
  using ModelType = Model<Matrix, SymmGroup>;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using ValueType = typename MPSType::value_type;
  using InitializerType = mps_initializer<Matrix, SymmGroup>;
  using LinearSystemSSSimulationType = SweepBasedLinearSystem<cmatrix, SymmGroup, StorageType, SweepOptimizationType::SingleSite>;
  using LinearSystemTSSimulationType = SweepBasedLinearSystem<cmatrix, SymmGroup, StorageType, SweepOptimizationType::TwoSite>;
  using PointerToSSSimulatorType = std::unique_ptr<LinearSystemSSSimulationType>;
  using PointerToTSSimulatorType = std::unique_ptr<LinearSystemTSSimulationType>;
  using ResultContainerType = std::map<std::pair<int, int>, MPSType>;

public:

  /** @brief Class constructor */
  FEASTSimulator(BaseParameters& parms, const ModelType& model, const LatticeType& lattice) 
    : currentIter(0), isSingleSite(true), parameters(parms) {
    // Retrieve simulation parameters
    numStates = parameters["feast_num_states"].as<int>();
    maxFeastIter = parameters["feast_max_iter"].as<int>();
    eMin = parameters["feast_emin"].as<double>();
    eMax = parameters["feast_emax"].as<double>();
    mMax = parameters["max_bond_dimension"].as<int>();
    feastThreshold = parameters["feast_convergence_threshold"].as<double>();
    numQuadraturePoint = parameters["feast_num_points"].as<int>();
    intModality = parameters["feast_integral_type"].as<std::string>();
    truncModality = parameters["feast_truncation_type"].as<std::string>();
    truncateEach = (truncModality == "each");
    initType = parameters["feast_init_type"].as<std::string>();
    printHeader();
    // Checks consistency of the input
    if (intModality != "half" && intModality != "full")
      throw std::runtime_error("Parameter [feast_integral_type] not recognized");
    if (truncModality != "each" && truncModality != "end")
      throw std::runtime_error("Parameter [feast_truncation_type] not recognized");
    // Generates the initial guess for the MPSs
    generateSeed(parms);
    initializeGuess(parms, model, lattice);
    quadPoints = FeastHelper::getQuadraturePoints(numQuadraturePoint);
    if (parameters["optimization"] == "twosite")
      isSingleSite = false;
  }

  /** @brief FEAST simulation (which is composed by multiple FEAST iterations) */
  void runFEAST(const MPOType& mpo) {

  }

  /** @brief Runs a single iteration of DMRG[FEAST] */
  void runFeastSimulation(const MPOType& mpo) {
    maquis::cout << " +====================" << std::endl;
    maquis::cout << "   FEAST iteration " <<  currentIter << std::endl;
    maquis::cout << " +====================" << std::endl;
    // Variable initialization
    std::vector<ComplexType> complexWeights;
    auto r = (eMax - eMin)/2.;
    auto r0 = (eMax + eMin - 2.*mpo.getCoreEnergy())/2.;
    // For each FEAST iteration, we have a loop over the number of quadrature points 
    // AND of the number of target states.
    for (int quadPoint = 0; quadPoint < numQuadraturePoint; quadPoint++) {
      double node = quadPoints[quadPoint].quadratureNode;
      double weight = quadPoints[quadPoint].quadratureWeight;
      double theta = (1.0 - node) * M_PI;
      if (intModality == "half")
        theta /= 2.;
      ComplexType rimag = r * std::exp(imagUnity * theta);
      ComplexType t = r0 + rimag;
      complexWeights.push_back(weight * rimag / 4.);
      maquis::cout << std::endl;
      maquis::cout << " == NEW QUADRATURE POINT ==" << std::endl;
      maquis::cout << std::endl;
      maquis::cout << " - Node: " << t << std::endl;
      maquis::cout << " - Weight: " << complexWeights[quadPoint] << std::endl;
      maquis::cout << std::endl;
      // Here there is a bit of code repetition because the pointer type is different for SS and TS
      maquis::cout << std::endl;
      if (isSingleSite) {
        for (int iGuess = 0; iGuess < numStates; iGuess++) {
          maquis::cout << std::endl;
          maquis::cout << " == Solving linear system for the guess " << iGuess << " ==" << std::endl;
          maquis::cout << " - Using seed: " << seedForInit[iGuess] << std::endl;
          auto mpsTmp = mpsGuess[iGuess];
          auto ssSimulator = std::make_unique<LinearSystemSSSimulationType>(mpsTmp, mpo, parameters, 0);
          ssSimulator->setShift(t);
          ssSimulator->runSweepSimulation();
          resultContainer.insert(std::make_pair(std::make_pair(iGuess, quadPoint), mpsTmp));
        }
      }
      else {
        for (int iGuess = 0; iGuess < numStates; iGuess++) {
          maquis::cout << " == Solving linear system for the guess " << iGuess << " ==" << std::endl;
          maquis::cout << " - Using seed: " << seedForInit[iGuess] << std::endl;
          auto mpsTmp = mpsGuess[iGuess];
          auto tsSimulator = std::make_unique<LinearSystemTSSimulationType>(mpsTmp, mpo, parameters, 0);
          tsSimulator->setShift(t);
          tsSimulator->runSweepSimulation();
          resultContainer.insert(std::make_pair(std::make_pair(iGuess, quadPoint), mpsTmp));
        }
      }
      maquis::cout << std::endl;
    }
    // Diagonalizes the Hamiltonian matrix in the FEAST subspace
    typename FeastHelper::FEASTPostProcessor<SymmGroup> postProcessor(resultContainer, numStates, numQuadraturePoint, complexWeights);
    postProcessor.solveEigenvalueProblem(mpo);
    resultContainer = postProcessor.performBackTransformation(mpo, mMax, truncateEach);
    postProcessor.printResults();
    energies = postProcessor.getEnergies();
    // Final update of the iteration counter
    currentIter += 1;
  }

  /** @brief Getter for the MPS guesses */
  auto getCurrentGuess(int iState) {
    return mpsGuess[iState];
  }

  /** @brief Getter for the quadrature points */
  auto getQuadraturePoints() const {
    return quadPoints;
  }

  /** @brief Getter for the FEAST energy */
  auto getEnergy(int iState) const {
    return energies[iState];
  }

private:

  /** @brief Generates the guess for FEAST */
  void initializeGuess(BaseParameters& parms, const ModelType& model, const LatticeType& lattice) {
    std::vector<std::string> initStates;
    bool needToWriteONV = (initType == "basis_state_generic" || initType == "basis_state_generic_const" || initType == "basis_state_generic_default");
    if (needToWriteONV) {
      std::string states = parms["feast_init_onv"].as<std::string>();
      initStates.resize(numStates);
      boost::split(initStates, states, boost::is_any_of("|"));
    }
    // Generates the guess MPS
    for (int iState = 0; iState < numStates; iState++) {
      auto parametersTmp = parms;
      parametersTmp["init_state"] = initType;
      parametersTmp["seed"] = seedForInit[iState];
      if (needToWriteONV)
        parametersTmp["init_basis_state"] = initStates[iState];
      mpsGuess.push_back(MPSType(lattice.size(), *(model.initializer(lattice, parametersTmp))));
    }
  }

  /** @brief Generates the seed for the random number generator */
  void generateSeed(BaseParameters& parms) {
    // Sets the seeds (to be used later)
    std::srand(parms["seed"]);
    seedForInit.resize(numStates);
    for (int iState = 0; iState < numStates; iState++)
      seedForInit[iState] = std::rand();
  }

  /** @brief Prints an header with the FEAST-specific parameters */
  void printHeader() const {
    maquis::cout << std::endl;
    maquis::cout << " =========================" << std::endl;
    maquis::cout << "  DMRGT[FEAST] SIMULATION " << std::endl;
    maquis::cout << " =========================" << std::endl;
    maquis::cout << std::endl;
    maquis::cout << " SIMULATION PARAMETERS: " << std::endl;
    maquis::cout << " - Number of targeted states: " << numStates << std::endl;
    maquis::cout << " - Lower bound for the complex contour integration: " << eMin << std::endl;
    maquis::cout << " - Uppwer bound for the complex contour integration: " << eMax << std::endl;
    maquis::cout << " - Number of quadrature points: " << numQuadraturePoint << std::endl;
    maquis::cout << " - Convergence threshold for FEAST: " << feastThreshold << std::endl;
    maquis::cout << " - Integration type: " << intModality << std::endl;
    maquis::cout << " - Truncation modality: " << truncModality << std::endl;
    maquis::cout << " - MPS guess type: " << initType << std::endl;
    maquis::cout << " - DMRG solver: " << ((isSingleSite) ? "single site" : "two site") << std::endl;
  }

  // -- Class members --
  BaseParameters parameters;                                     // Parameter container
  int currentIter;                                               // Index of the current FEAST iteration.
  int numStates;                                                 // Number of states to be targeted.
  int maxFeastIter;                                              // Maximum number of FEAST iterations.
  int mMax;                                                      // Maximum value of the bond dimension.
  double eMin, eMax;                                             // Lower and upper bound for the complex contour integral.
  double feastThreshold;                                         // Threshold to assess the convergence of DMRG[FEAST].
  int numQuadraturePoint;                                        // Number of quadrature point.
  std::string intModality;                                       // "Full" for the full circle integration, "half" for the half-circle one.
  std::string truncModality;                                     // "Each" if the MPS must be truncated after each sum, "end" if the truncation must be done only at the end.
  std::string initType;                                          // Initialization strategy for each guess.
  std::vector<MPSType> mpsGuess;                                 // Stores the current guess for hte FEAST procedure.
  std::vector<int> seedForInit;                                  // Seed for random initialization.
  std::vector<typename FeastHelper::QuadraturePoint> quadPoints; // Vector with the quadrature points and weight.
  bool isSingleSite;                                             // If true, runs a single-site calculation, otherwise runs a two-sites one.
  bool truncateEach;                                             // If true, truncates the MPS after each sum.
  ResultContainerType resultContainer;                           // Member that stores the result of each linear system.
  std::vector<double> energies;                                  // FEAST energies at the current iteration.
  // Constexpr for the imaginary unit
  static constexpr ComplexType imagUnity = ComplexType(0., 1.);
};

template<class SymmGroup>
constexpr typename FEASTSimulator<SymmGroup>::ComplexType FEASTSimulator<SymmGroup>::imagUnity;

#endif // FEAST_SIMULATOR
