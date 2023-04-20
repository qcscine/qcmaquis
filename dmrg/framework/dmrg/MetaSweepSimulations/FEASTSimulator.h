/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef FEAST_SIMULATOR
#define FEAST_SIMULATOR

#include <omp.h>
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
  using PostProcessorType = typename FeastHelper::FEASTPostProcessor<SymmGroup>;

public:

  /** @brief Class constructor */
  FEASTSimulator(BaseParameters& parms, const ModelType& model, const LatticeType& inputLattice, const MPOType& mpo)
    : currentIter(0), isSingleSite(true), parameters(parms), mpo_(mpo), lattice(inputLattice), calculateExactError(false),
      model_(model), feastMPSs(), calculateVariance(false)
  {
    // Retrieve simulation parameters
    verbose_ = (parameters["feast_verbose"] == "yes");
    printTimings_ = (parameters["feast_print_timings"] == "yes");
    numStates = parameters["feast_num_states"].as<int>();
    maxFeastIter = parameters["feast_max_iter"].as<int>();
    eMin = parameters["feast_emin"].as<double>();
    eMax = parameters["feast_emax"].as<double>();
    mMax = parameters["max_bond_dimension"].as<int>();
    feastThresholdEnergy = parameters["feast_energy_convergence_threshold"].as<double>();
    feastThresholdOverlap = parameters["feast_overlap_convergence_threshold"].as<double>();
    numQuadraturePoint = parameters["feast_num_points"].as<int>();
    intModality = parameters["feast_integral_type"].as<std::string>();
    truncModality = parameters["feast_truncation_type"].as<std::string>();
    truncateEach = (truncModality == "each");
    initType = parameters["init_type"].as<std::string>();
    if (parameters["linsystem_exact_error"] == "yes")
      calculateExactError = true;
    calculateVariance = (parameters["feast_calculate_standard_deviation"] == "yes");
    // Checks consistency of the input
    if (intModality != "half" && intModality != "full")
      throw std::runtime_error("Parameter [feast_integral_type] not recognized");
    if (truncModality != "each" && truncModality != "end")
      throw std::runtime_error("Parameter [feast_truncation_type] not recognized");
    // Generates the initial guess for the MPSs
    generateSeed(parms);
    initializeGuess(parms, model_);
    quadPoints = FeastHelper::getQuadraturePoints(numQuadraturePoint);
    this->generateComplexQuadrature();
    if (parameters["optimization"] == "twosite")
      isSingleSite = false;
    postProcessor = std::make_unique<PostProcessorType>(numStates, numQuadraturePoint, complexWeights, model_, lattice,
                                                        parameters, eMin, eMax);
    resultContainer = std::make_shared<ResultContainerType>();
    printHeader();
  }

  /** @brief FEAST simulation (which is composed by multiple FEAST iterations) */
  void runFEAST() {
    maquis::cout << std::endl;
    maquis::cout << " +===============================+" << std::endl;
    maquis::cout << "   STARTING NEW FEAST SIMULATION" << std::endl;
    maquis::cout << " +===============================+" << std::endl;
    maquis::cout << std::endl;
    bool converged=false;
    while (!converged) {
      maquis::cout << " ---------------------" << std::endl;
      maquis::cout << "   FEAST iteration " <<  currentIter << std::endl;
      maquis::cout << " ---------------------" << std::endl;
      this->runFeastSimulation();
      // Prints information about convergence
      if (currentIter > 1) {
        // Checks overlap convergence criterion
        auto overlapVector = this->getOverlapMatrix();
        double overlapDifference = std::accumulate(overlapVector.begin(), overlapVector.end(), 0., [](int sum, double overlap) {
          return sum + (1.-std::abs(overlap));
        });
        // Checks the energy convergence criterion
        auto energyDifference = postProcessor->getOverallEnergyVariation();
        // Final printing
        maquis::cout << std::endl;
        maquis::cout << " +-- CONVERGENCE CHECK --+" << std::endl;
        maquis::cout << std::endl;
        maquis::cout << "  - Iteration: " << currentIter << std::endl;
        maquis::cout << std::scientific << std::setprecision(16);
        maquis::cout << "  - Overall MPS variation: " << overlapDifference << std::endl;
        maquis::cout << "  - Overall energy variation: " << energyDifference << std::endl;
        bool hasReachedConvergence = ((overlapDifference < feastThresholdOverlap) && (energyDifference < feastThresholdEnergy));
        converged = (currentIter == maxFeastIter) || hasReachedConvergence;
        maquis::cout << std::endl;
        if (converged) {
          if (currentIter == maxFeastIter && !hasReachedConvergence) {
            maquis::cout << " --> MAXIMUM NUMBER OF FEAST ITERATION REACHED, STOPPING DMRG[FEAST]" << std::endl;
          }
          else {
            maquis::cout << " --> CONVERGENCE REACHED" << std::endl;
          }
        }
        else {
          maquis::cout << " --> CONVERGENCE NOT REACHED, WILL START NEW ITERATION" << std::endl;
        }
        maquis::cout << std::endl;
      }
      else {
        converged = (currentIter == maxFeastIter);
      }
    }
  }

  /** @breif Getter for the guess */
  auto getCurrentGuess(int iState) const {
    if (iState >= mpsGuess.size())
      throw std::runtime_error("FEAST: guess index not available");
    else
      return mpsGuess[iState];
  }

  /** @brief Getter for the MPS guesses */
  auto getCurrentEigenvalues() {
    return screenedMPSs;
  }

  /** @brief Getter for the MPS guesses */
  auto getCurrentEigenvalue(int iState) {
    return screenedMPSs->operator[](iState);
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

  /** @brief Runs a single iteration of DMRG[FEAST] */
  void runFeastSimulation() {
    // Variable initialization
    // Note that, if this is not the first FEAST iteration, then the feastMPSs vector is defined
    // and we can generate the initial guesses from the results of the previous FEAST iteration. 
    if (feastMPSs)
      for (int iState = 0; iState < feastMPSs->size(); iState++)
        mpsGuess[iState] = feastMPSs->operator[](iState);
    // For each FEAST iteration, we have a loop over the number of quadrature points
    // *and* of the number of target states.
    auto initialTime = std::chrono::high_resolution_clock::now();
    //
#pragma omp parallel for collapse(2)
    for (int quadPoint = 0; quadPoint < numQuadraturePoint; quadPoint++) {
      for (int iGuess = 0; iGuess < numStates; iGuess++) {
        auto localParameters = parameters;
        auto mpsTmp = mpsGuess[iGuess];
        // Here there is a bit of code repetition because the pointer type is different for SS and TS.
        if (isSingleSite) {
          auto ssSimulator = std::make_unique<LinearSystemSSSimulationType>(mpsTmp, mpo_, localParameters, model_, lattice, verbose_);
          ssSimulator->setShift(complexNodes[quadPoint]);
          auto initialInnerTime = std::chrono::high_resolution_clock::now();
          ssSimulator->runSweepSimulation();
          auto finalInnerTime = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration<double>(finalInnerTime - initialInnerTime)).count();
#pragma omp critical (PrintResults)
          {
            if (!verbose_) {
              printLinearSystemHeader(complexNodes[quadPoint], complexWeights[quadPoint], iGuess, duration);
              ssSimulator->printSummary();
            }
          }
        }
        else {
          auto tsSimulator = std::make_unique<LinearSystemTSSimulationType>(mpsTmp, mpo_, localParameters, model_, lattice, verbose_);
          tsSimulator->setShift(complexNodes[quadPoint]);
          auto initialInnerTime = std::chrono::high_resolution_clock::now();
          tsSimulator->runSweepSimulation();
          auto finalInnerTime = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration<double>(finalInnerTime - initialInnerTime)).count();
#pragma omp critical (PrintResults)
          {
            if (!verbose_) {
              printLinearSystemHeader(complexNodes[quadPoint], complexWeights[quadPoint], iGuess, duration);
              tsSimulator->printSummary();
            }
          }
        }
#pragma omp critical (UpdateOfResults)
        {
          // Final update of the results
          resultContainer->insert(std::make_pair(std::make_pair(iGuess, quadPoint), mpsTmp));
        }
      }
    }
    auto finalTime = std::chrono::high_resolution_clock::now();
    if (printTimings_) {
      double elapsed = std::chrono::duration<double>(finalTime - initialTime).count();
      maquis::cout << " --> Time elapsed to solve the FEAST problem " << elapsed << "." << std::endl;
    }
    // Diagonalizes the Hamiltonian matrix in the FEAST subspace
    postProcessor->updateContainer(resultContainer);
    postProcessor->solveEigenvalueProblem(mpo_);
    feastMPSs = postProcessor->performBackTransformation(mpo_, mMax, truncateEach);
    if (screenedMPSs)
      screenedMPSsPrevious = screenedMPSs;
    screenedMPSs = postProcessor->getScreenedMPSs();
    postProcessor->printResults();
    energies = postProcessor->getScreenedEnergies();
    // Final update of the iteration counter
    currentIter += 1;
  }

  /** @brief Get overlap variation */
  auto getOverlapMatrix() const {
    auto overlaps = std::vector<double>(screenedMPSs->size(), 0);
    maquis::cout << std::endl;
    maquis::cout << " == OVERLAP CONVERGENCE CHECK == " << std::endl;
    maquis::cout << std::endl;
    // Finds, among the previous MPSs, the ones which have the largest overlap
    // with the current FEAST MPSs.
    int iMPS = 0;
    for (const auto& iCurrent: *screenedMPSs) {
      int iFound = 0, iRef = 0;
      for (const auto& iPrevious: *screenedMPSsPrevious) {
        auto localOverlap = std::abs(overlap(iCurrent, iPrevious))/std::sqrt(norm(iCurrent)*norm(iPrevious));
        if (localOverlap > overlaps[iMPS]) {
          overlaps[iMPS] = localOverlap;
          iRef = iFound;
        }
        iFound += 1;
      }
      maquis::cout << " State " << iMPS << " has the largest overlap with state " << iRef << " of the previous iteration, overlap = " 
                   << overlaps[iMPS] << std::endl;
      iMPS += 1;
    }
    maquis::cout << std::endl;
    return overlaps;
  }

  /** @brief Generates the complex-valued nodes and weights */
  void generateComplexQuadrature() {
    auto r = (eMax - eMin)/2.;
    auto r0 = (eMax + eMin - 2.*mpo_.getCoreEnergy())/2.;
    for (int quadPoint = 0; quadPoint < numQuadraturePoint; quadPoint++) {
      double node = quadPoints[quadPoint].quadratureNode;
      double weight = quadPoints[quadPoint].quadratureWeight;
      double theta = (1.0 - node) * M_PI;
      if (intModality == "half")
        theta /= 2.;
      ComplexType rimag = r * std::exp(imagUnity * theta);
      complexNodes.push_back(r0 + rimag);
      complexWeights.push_back(weight * rimag / 4.);
    }
  }

  /** @brief Generates the guess for FEAST */
  void initializeGuess(BaseParameters& parms, const ModelType& model) {
    bool needToWriteONV = (initType == "basis_state_generic" || initType == "hf");
    std::vector<std::string> specifiedStates;
    int numSpecifiedStates = 0;
    if (needToWriteONV) {
      std::string states = parms["init_basis_state"].as<std::string>();
      boost::split(specifiedStates, states, boost::is_any_of("|"));
      numSpecifiedStates = specifiedStates.size();
      if (numSpecifiedStates < 1 || numSpecifiedStates > numStates){
        throw std::runtime_error("You should specify at least one and at most num_states init_onv's if init_type is set to basis_state_generic/hf");
      }
      if (numSpecifiedStates != numStates) {
        maquis::cout << "WARNING! Not all feast states have been provided an ONV for initialization, so the remaining ones will be initialized with generic_default" << std::endl;
      }
    }
    // Generates the guess MPS
    for (int iState = 0; iState < numStates; iState++) {
      auto parametersTmp = parms;
      parametersTmp.set("seed", seedForInit[iState]);
      parametersTmp.set("init_type", initType);
      if (needToWriteONV) {
        if (iState < numSpecifiedStates) {
          if (parametersTmp["MODEL"] == "quantum_chemistry")
            parametersTmp.set("hf_occ", specifiedStates[iState]);
          else
            parametersTmp.set("init_basis_state", specifiedStates[iState]); // initialize the specified states with the provided ONVs
        } else { // the rest of the states are not specified
          if (parametersTmp["MODEL"] == "quantum_chemistry")
            parametersTmp.set("init_type", "const "); // initialize the remaining electronic states with const
          else
            parametersTmp.set("init_type", "basis_state_generic_default"); // initialize the remaining vibrational states with generic_default
        }
      }
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
    maquis::cout << " - Maximum number of FEAST iterations: " << maxFeastIter << std::endl;
    maquis::cout << " - Number of targeted states: " << numStates << std::endl;
    maquis::cout << " - Lower bound for the complex contour integration: " << eMin << std::endl;
    maquis::cout << " - Uppwer bound for the complex contour integration: " << eMax << std::endl;
    maquis::cout << " - Number of quadrature points: " << numQuadraturePoint << std::endl;
    maquis::cout << " - Energy convergence threshold for FEAST: " << feastThresholdEnergy << std::endl;
    maquis::cout << " - Overlap convergence threshold for FEAST: " << feastThresholdOverlap << std::endl;
    maquis::cout << " - Integration type: " << intModality << std::endl;
    maquis::cout << " - Truncation modality: " << truncModality << std::endl;
    maquis::cout << " - MPS guess type: " << initType << std::endl;
    maquis::cout << " - DMRG solver: " << ((isSingleSite) ? "single site" : "two site") << std::endl;
    if (calculateVariance)
      maquis::cout << " - Calculating variance for each eigenpair." << std::endl;
  }

  /** @brief Prints the header for the solution of a given linear system */
  static void printLinearSystemHeader(ComplexType node, ComplexType weight, int iGuess, int time) {
    maquis::cout << std::endl;
    maquis::cout << " == NEW FEAST LINEAR SYSTEM ==" << std::endl;
    maquis::cout << std::endl;
    maquis::cout << " - Node: " << node << std::endl;
    maquis::cout << " - Weight: " << weight << std::endl;
    maquis::cout << " - Guess number: " << iGuess << std::endl;
    maquis::cout << " - Time elapsed in linear system solver: " << time << " milliseconds" << std::endl;
    maquis::cout << std::endl;
  }

  // -- Class members --
  BaseParameters parameters;                                                           // Parameter container
  int currentIter;                                                                     // Index of the current FEAST iteration.
  int numStates;                                                                       // Number of states to be targeted.
  int maxFeastIter;                                                                    // Maximum number of FEAST iterations.
  int mMax;                                                                            // Maximum value of the bond dimension.
  double eMin, eMax;                                                                   // Lower and upper bound for the complex contour integral.
  double feastThresholdEnergy, feastThresholdOverlap;                                  // Threshold to assess the convergence of DMRG[FEAST] (both energy and overlap).
  int numQuadraturePoint;                                                              // Number of quadrature point.
  std::string intModality;                                                             // "Full" for the full circle integration, "half" for the half-circle one.
  std::string truncModality;                                                           // "Each" if the MPS must be truncated after each sum, "end" if the truncation must be done only at the end.
  std::string initType;                                                                // Initialization strategy for each guess.
  std::vector<MPSType> mpsGuess;                                                       // Stores the current guess for hte FEAST procedure.
  std::shared_ptr<std::vector<MPSType>> feastMPSs, screenedMPSs, screenedMPSsPrevious; // Final, back-transformed FEAST MPSs, and their screened counterpart.
  std::vector<int> seedForInit;                                                        // Seed for random initialization.
  std::vector<typename FeastHelper::QuadraturePoint> quadPoints;                       // Vector with the quadrature points and weight.
  bool isSingleSite;                                                                   // If true, runs a single-site calculation, otherwise runs a two-sites one.
  bool truncateEach;                                                                   // If true, truncates the MPS after each sum.
  bool calculateExactError;                                                            // If true, calculates the exact error associated with the linear system.
  bool calculateVariance;                                                              // If true, calculates the variance at the end of each FEAST iteration.
  std::shared_ptr<ResultContainerType> resultContainer;                                // Member that stores the result of each linear system.
  std::vector<double> energies;                                                        // FEAST energies at the current iteration.
  std::unique_ptr<PostProcessorType> postProcessor;                                    // Class managing FEAST postprocessing.
  std::vector<ComplexType> complexNodes, complexWeights;                               // Quadrature rule for the complex circle.
  const MPOType& mpo_;                                                                 // Matrix product operator
  const LatticeType& lattice;                                                          // DMRG lattice object.
  const ModelType& model_;                                                             // Model object.
  bool verbose_, printTimings_;                                                        // Verbosity flags.
  // Constexpr for the imaginary unit
  static constexpr ComplexType imagUnity = ComplexType(0., 1.);
};

template<class SymmGroup>
constexpr typename FEASTSimulator<SymmGroup>::ComplexType FEASTSimulator<SymmGroup>::imagUnity;

#endif // FEAST_SIMULATOR