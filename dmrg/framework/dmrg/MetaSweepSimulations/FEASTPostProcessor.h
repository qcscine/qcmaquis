/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef FEAST_HELPER_CLASS
#define FEAST_HELPER_CLASS

#include <stdexcept>
#include <vector>
#include "alps/numeric/matrix/algorithms.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"

namespace FeastHelper {

/** @brief Enum class representing whether a root is accepted or not */
enum class EigenvalueSelection {Accepted, NotInInterval, HighVariance };

/** @brief Class devoted to the post-processing of the FEAST data */
template <class SymmGroup>
class FEASTPostProcessor {
public:
  using ChargeType = typename SymmGroup::charge;
  using LatticeType = Lattice;
  using Matrix = cmatrix;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using ModelType = Model<Matrix, SymmGroup>;
  using RealValueType = double;
  using RealMatrixType = alps::numeric::matrix<RealValueType>;
  using RealVectorType = alps::numeric::vector<RealValueType>;
  using ComplexNumber = typename std::complex<RealValueType>;
  using ComplexMatrixType = alps::numeric::matrix<ComplexNumber>;
  using DiagonalMatrixType = typename alps::numeric::associated_real_diagonal_matrix<ComplexMatrixType>::type;
  using ComplexVectorType = alps::numeric::vector<ComplexNumber>;
  using ResultContainerType = std::map<std::pair<int, int>, MPSType>;
  using VectorOfMPSs = std::vector<MPSType>;
  using MatrixOfMPSs = std::map<std::pair<int, int>, MPSType >;

  FEASTPostProcessor(int numberOfStates, int numberOfQuadrature, const std::vector<ComplexNumber>& w, const ModelType& inputModel,
                     const LatticeType& inputLattice, BaseParameters& parms_, double eMin_, double eMax_)
    : nStates(numberOfStates), nQuad(numberOfQuadrature), weights(w), model(inputModel), lattice(inputLattice), calculateStandardDeviation(false),
      eMin(eMin_), eMax(eMax_), parms(parms_), screenedEnergies()
  {
    energies = std::vector<double>(nStates, 0);
    energiesPrev = std::vector<double>(nStates, 0);
    truncatedEnergy = std::vector<double>(nStates, 0);
    standardDeviations = std::vector<double>(nStates, 0);
    accepted = std::vector<EigenvalueSelection>(nStates, EigenvalueSelection::Accepted);
    totalQN = model.total_quantum_numbers(parms);
    calculateStandardDeviation = (parms["feast_calculate_standard_deviation"] == "yes");
    screenStandardDeviation = parms.is_set("feast_standard_deviation_threshold");
    if (screenStandardDeviation)
      standardDeviationScreening = parms["feast_standard_deviation_threshold"].as<double>();
  }

  /** @brief Updates teh mps container */
  void updateContainer(std::shared_ptr<ResultContainerType> container) {
    mpsContainer = container;
  }

  /**
   * @brief Diagonalizes the Hamiltonian in the FEAST subspace.
   * @return std::vector<double> eigenvalues of the FEAST Hamiltonian.
   */
  void solveEigenvalueProblem(const MPOType& mpo) {
    // -- Hamiltonian matrix construction --
    // auto Hvec = std::vector<cmat_type>(omp_get_max_threads(), cmat_type::Zero(n_states, n_states));
    // auto Bvec = std::vector<cmat_type>(omp_get_max_threads(), cmat_type::Zero(n_states, n_states));
    ComplexMatrixType H = ComplexMatrixType(nStates, nStates, 0.);
    ComplexMatrixType B = ComplexMatrixType(nStates, nStates, 0.);
    maquis::cout << std::endl;
    maquis::cout << " +------------------------------+" << std::endl;
    maquis::cout << "  FEAST SUBSPACE DIAGONALIZATION" << std::endl;
    maquis::cout << " +------------------------------+" << std::endl;
    maquis::cout << std::endl;
    // #pragma omp parallel for collapse(2)
    // Note that here we do not first sum the MPS and then calculate the expectation value, because this
    // would lead to a very large MPS. We instead sum the expectation values directly.
    for (int i = 0; i < nStates; i++) {
      for (int iQ = 0; iQ < nQuad; iQ++) {
        auto& mpsCopyI = mpsContainer->operator[](std::make_pair(i, iQ));
        for (int j = 0; j < nStates; j++) {
          for (int jQ = 0; jQ < nQuad; jQ++) {
            auto& mpsCopyJ = mpsContainer->operator[](std::make_pair(j, jQ));
            H(i, j) += expval(mpsCopyI, mpsCopyJ, mpo)*std::conj(weights[iQ])*weights[jQ];
            B(i, j) += overlap(mpsCopyI, mpsCopyJ)*std::conj(weights[iQ])*weights[jQ];
          }
        }
      }
    }
    // maquis::cout << " Hamiltonian matrix in the FEAST subspace" << std::endl;
    // maquis::cout << H << std::endl;
    // maquis::cout << " Overlap matrix of the FEAST subspace" << std::endl;
    // maquis::cout << B << std::endl;
    // for (int i = 0; i < nStates; i++)
    //   normVector[i] = std::sqrt(std::real(B(i, i)));
    // for (int i = 0; i < nStates; i++) {
    //   for (int j = 0; j < nStates; j++) {
    //     H(i, j) /= normVector[i]*normVector[j];
    //     B(i, j) /= normVector[i]*normVector[j];
    //   }
    // }
    // == Matrix diagonalization ==
    auto zeroComplex = ComplexNumber(0., 0.);
    ComplexMatrixType leftEigenVectors(nStates, nStates, zeroComplex), rightEigenVectors(nStates, nStates, zeroComplex);
    ComplexVectorType alphaVec(nStates, 0.), betaVec(nStates, zeroComplex);
    alps::numeric::ggev(H, B, alphaVec, betaVec, leftEigenVectors, rightEigenVectors, thresholdForRank_);
    // // Calculates the overlap in the right eigenvectors basis
    // ComplexMatrixType overlapEigenvectorBasis(nStates, nStates);
    // gemm(adjoint(rightEigenVectors), rightEigenVectors, overlapEigenvectorBasis);
    // RealVectorType overlapEigenValues(nStates);
    // heev(overlapEigenvectorBasis, overlapEigenValues);
    // std::cout << overlapEigenValues << std::endl;
    // Retrieves energies
    energiesPrev = energies;
    screenedEnergiesPrev = screenedEnergies;
    for (int iState = 0; iState < nStates; iState++) {
      if (std::imag(alphaVec[iState]) > 1.0E-10 || std::imag(betaVec[iState]) > 1.0E-10)
        maquis::cout << " WARNING: Energy of the " << iState << "-th state has a non-negligible imaginary part" << std::endl;
      if (std::abs(betaVec[iState]) > thresholdForRank_) {
        energies[iState] = maquis::real(alphaVec[iState]/betaVec[iState]);
        accepted[iState] = (energies[iState] > eMin && energies[iState] < eMax) ? EigenvalueSelection::Accepted
                                                                                : EigenvalueSelection::NotInInterval;
      }
      else {
        energies[iState] = 0.;
        accepted[iState] = EigenvalueSelection::NotInInterval;
      }
    }
    // Final copy of the eigenvectors
    feastEigenVectors = rightEigenVectors;
  };

  /**
   * @brief Back-transformation of the MPS
   *
   * The eigenvectors of the FEAST Hamiltonian are used to back-transform the MPS
   * and provide new guess for FEAST iteration.
   *
   * @param mMax maximum bond dimension
   * @param truncEach if true, truncates after each sum between MPSs.
   */
  auto performBackTransformation(const MPOType& mpo, int mMax, bool truncEach) {
    // Generates the MPS files for the new FEAST iteration
    // using MatrixOfMPSs = Eigen::Matrix< MPS<cMatrix, SymmGroup>, -1, -1>;
    // MatrixOfMPSs mps_transf(n_states, n_states);
    MatrixOfMPSs mpsTransformed;
    // Variable definition
    auto refNorm = ietl::two_norm(mpsContainer->begin()->second[0]);
    currentFEASTMPSs = std::make_shared<VectorOfMPSs>(nStates);
    for (auto& iMPS: *mpsContainer)
      iMPS.second[0] /= refNorm;
    // Actual back-transformation
    //#pragma omp parallel for collapse(2)
    for (int iOutput = 0; iOutput < nStates; iOutput++) {
      if (accepted[iOutput] == EigenvalueSelection::Accepted) {
        for (int iInput = 0; iInput < nStates; iInput++) {
          // maquis::cout << "(" << iInput << "," << iOutput << ") = " << feastEigenVectors(iInput, iOutput) << std::endl;
          if (std::abs(feastEigenVectors(iInput, iOutput)) > thresholdForRank_) {
            for (int iQuad = 0; iQuad < nQuad; iQuad++) {
              // maquis::cout << " - iQuad = " << weights[iQuad] << std::endl;
              MPSType mpsToAdd = mpsContainer->operator[](std::make_pair(iInput, iQuad));
              auto scalingFactor = feastEigenVectors(iInput, iOutput)*weights[iQuad]; // /normVector[iInput];
              // std::cout << scalingFactor << std::endl;
              mpsToAdd.scaleByScalar(scalingFactor);
              if (iQuad == 0) {
                mpsTransformed[std::make_pair(iOutput, iInput)] = mpsToAdd;
              }
              else {
                if (std::abs(scalingFactor) > thresholdForRank_) {
                  if (truncEach)
                    mpsTransformed[std::make_pair(iOutput, iInput)] = joinAndTruncate(mpsTransformed[std::make_pair(iOutput, iInput)], mpsToAdd, mMax);
                  else
                    mpsTransformed[std::make_pair(iOutput, iInput)] = join(mpsTransformed[std::make_pair(iOutput, iInput)], mpsToAdd);
                }
              }
            }
            if (!truncEach)
              mpsTransformed[std::make_pair(iOutput, iInput)] = compression::l2r_compress(mpsTransformed[std::make_pair(iOutput, iInput)], mMax, thresholdForRank_, false);
          }
        }
      }
    }
    //#pragma omp parallel for
    for (int iOutput = 0; iOutput < nStates; iOutput++) {
      // If the state has been accepted,
      if (accepted[iOutput] == EigenvalueSelection::Accepted) {
        // Finds the first non-zero MPSs
        bool found=false;
        int iFirstInput = 0;
        while (!found) {
          auto position = mpsTransformed.find(std::make_pair(iOutput, iFirstInput));
          if (position != mpsTransformed.end())
            found = true;
          else
            iFirstInput += 1;
        }
        currentFEASTMPSs->operator[](iOutput) = mpsTransformed[std::make_pair(iOutput, iFirstInput)];
        for (int iInput = iFirstInput+1; iInput < nStates; iInput++) {
          auto key = std::make_pair(iOutput, iInput);
          if (mpsTransformed.find(key) != mpsTransformed.end()) {
            if (truncEach)
              currentFEASTMPSs->operator[](iOutput) = joinAndTruncate(currentFEASTMPSs->operator[](iOutput), mpsTransformed[key], mMax);
            else
              currentFEASTMPSs->operator[](iOutput) = join(currentFEASTMPSs->operator[](iOutput), mpsTransformed[key]);
          }
        }
        currentFEASTMPSs->operator[](iOutput).normalize_right();
        if (!truncEach)
          currentFEASTMPSs->operator[](iOutput) = compression::l2r_compress(currentFEASTMPSs->operator[](iOutput), mMax, 1.0E-16);
        truncatedEnergy[iOutput] = maquis::real(expval(currentFEASTMPSs->operator[](iOutput), mpo)/norm(currentFEASTMPSs->operator[](iOutput)));
        // If requested, calculates the standard deviations
        if (calculateStandardDeviation)
          standardDeviations[iOutput] = this->getStandardDeviation(mpo, currentFEASTMPSs->operator[](iOutput), mMax);
      }
      // If the state has not been accepted, just regenerates a random vector
      else {
        currentFEASTMPSs->operator[](iOutput) = generateRandomMPS();
        currentFEASTMPSs->operator[](iOutput).normalize_right();
      }
    }
    // If a variance-based screening has been requested, overrides again with a random value
    if (screenStandardDeviation && calculateStandardDeviation) {
      for (int iState = 0; iState < nStates; iState++) {
        if (std::abs(standardDeviations[iState]) > standardDeviationScreening && accepted[iState] == EigenvalueSelection::Accepted) {
          accepted[iState] = EigenvalueSelection::HighVariance;
          currentFEASTMPSs->operator[](iState) = generateRandomMPS();
          currentFEASTMPSs->operator[](iState).normalize_right();
        }
      }
    }
    // Finalizes
    screenedEnergies.clear();
    for (int iState = 0; iState < energies.size(); iState++)
      if (accepted[iState] == EigenvalueSelection::Accepted)
        screenedEnergies.push_back(energies[iState]);
    return currentFEASTMPSs;
  }

  /** @brief Prints the results of the FEAST calculation */
  void printResults() const {
    maquis::cout << " +----------------------------------------------------------+----------------------+" << std::endl;
    maquis::cout << " |   State    |      Old energy      |      New energy      |   Truncated energy   |" << std::endl;
    maquis::cout << " +----------------------------------------------------------+----------------------+" << std::endl;
    for (int iState = 0; iState < energies.size(); iState++) {
      switch (accepted[iState]) {
        case EigenvalueSelection::Accepted:
          maquis::cout << std::setw(13) << std::internal << iState
                       << std::setw(23) << std::right << std::fixed << std::setprecision(8) << energiesPrev[iState]
                       << std::setw(23) << std::right << std::fixed << std::setprecision(8) << energies[iState]
                       << std::setw(23) << std::right << std::fixed << std::setprecision(8) << truncatedEnergy[iState]
                       << "  --> ROOT ACCEPTED";
          break;
        case EigenvalueSelection::NotInInterval:
          maquis::cout << std::setw(13) << std::internal << iState
                       << std::setw(23) << std::right << std::fixed << std::setprecision(8) << energiesPrev[iState]
                       << std::setw(23) << std::right << std::fixed << std::setprecision(8) << energies[iState]
                       << "          #########      --> ROOT NOT ACCEPTED: energy outside boundary";
          break;
        case EigenvalueSelection::HighVariance:
          maquis::cout << std::setw(13) << std::internal << iState
                       << std::setw(23) << std::right << std::fixed << std::setprecision(8) << energiesPrev[iState]
                       << std::setw(23) << std::right << std::fixed << std::setprecision(8) << energies[iState]
                       << "          #########      --> ROOT NOT ACCEPTED: variance too high";
          break;
      }
      maquis::cout << std::endl;
    }
    maquis::cout << " +---------------------------------------------------------------------------------+" << std::endl;
    maquis::cout << std::endl;
    // If requested, prints also the standard deviations
    if (calculateStandardDeviation) {
      maquis::cout << " +---------------------------------------------+" << std::endl;
      maquis::cout << " |   State    |    Energy standard deviation   |" << std::endl;
      maquis::cout << " +---------------------------------------------+" << std::endl;
      for (int iState = 0; iState < energies.size(); iState++) {
        if (accepted[iState] == EigenvalueSelection::Accepted) {
          maquis::cout << std::setw(13) << std::internal << iState
                       << std::setw(23) << std::right << std::fixed << std::setprecision(8) << standardDeviations[iState]
                       << std::endl;
        } else {
          maquis::cout << std::setw(13) << std::internal << iState
            << "             #########             --> ROOT NOT ACCEPTED: variance not calculated" << std::endl;
        }
      }
      maquis::cout << " +---------------------------------------------+" << std::endl;
      maquis::cout << std::endl;
    }
  }

  /** @brief Getter for the vibrational energy */
  auto getEnergies() const {
    return energies;
  }

  /** @brief Getter for the vibrational energy */
  auto getScreenedEnergies() const {
    return screenedEnergies;
  }

  /** @brief Getter for the screened MPSs */
  auto getScreenedMPSs() const {
    auto screenedMPS = std::make_shared<VectorOfMPSs>();
    for (int iState = 0; iState < currentFEASTMPSs->size(); iState++)
      if (accepted[iState] == EigenvalueSelection::Accepted)
        screenedMPS->push_back(currentFEASTMPSs->operator[](iState));
    return screenedMPS;
  }

  /** @brief Gets the overall energy variation */
  auto getOverallEnergyVariation() const {
    auto overallSum = std::accumulate(screenedEnergies.begin(), screenedEnergies.end(), 0.);
    auto oldOverallSum = std::accumulate(screenedEnergiesPrev.begin(), screenedEnergiesPrev.end(), 0.);
    return std::abs(overallSum-oldOverallSum);
  }

private:

  /** @brief Generates a random MPS */
  auto generateRandomMPS() {
    auto tmpParms = parms;
    tmpParms.set("init_type", "default");
    tmpParms.set("seed", std::rand());
    return MPSType(lattice.size(), *(model.initializer(lattice, tmpParms)));
  }

  /** @brief Static method to calculate determinant */
  static ComplexNumber calculateDeterminant(ComplexMatrixType inputMatrix) {
    int numRows = num_rows(inputMatrix);
    // Note that here we use long int for coherence with lapack
    std::vector<long int> ipiv(numRows);
    int info = boost::numeric::bindings::lapack::getrf(inputMatrix, ipiv);
    if (info != 0)
      throw std::runtime_error("Error in LU decomposition");
    ComplexNumber det = ComplexNumber(1., 0.);
    for (int iElement = 0; iElement < numRows; iElement++)
      det *= inputMatrix(iElement, iElement);
    return det;
  }

  /** @brief Standard deviation calculator -- needed to screen the FEAST eigenfunctions */
  auto getStandardDeviation(const MPOType& mpo, const MPSType& inputMPS, int mMax) {
    auto traitClass = MPOTimesMPSTraitClass<Matrix, SymmGroup>(inputMPS, model, lattice, totalQN, mMax);
    auto outputMPS = traitClass.applyMPO(mpo);
    auto squaredEnergy = (overlap(outputMPS, outputMPS) + 2.*maquis::real(overlap(outputMPS, inputMPS)*mpo.getCoreEnergy())
                           + overlap(inputMPS, inputMPS)*std::norm(mpo.getCoreEnergy()))/norm(inputMPS); // this is <H^2>
    auto energy = expval(inputMPS, mpo)/norm(inputMPS);
    auto energySquared = std::norm(energy); // this is <H>^2
    return std::sqrt(maquis::real(squaredEnergy - energySquared));
  }

  // -- Class members --
  BaseParameters& parms;                                         // Parameter container
  std::shared_ptr<ResultContainerType> mpsContainer;             // Data structure storing the result of the FEAST linear systems.
  int nStates, nQuad;                                            // FEAST-specific integer parameters.
  std::vector<double> energies, energiesPrev, screenedEnergies,  // Energies (possibly screened) at the current and previous iteration.
    screenedEnergiesPrev;
  std::vector<EigenvalueSelection> accepted;                     // Eigenpairs that are accepted.
  std::vector<double> truncatedEnergy, standardDeviations;       // FEAST-specific double parameters for checks.
  ComplexMatrixType feastEigenVectors;                           // FEAST --> eigenvalues transformation matrix.
  static constexpr int thresholdForRank_ = 1.0E-10;              // Threshold for rank.
  RealVectorType eigenValues;                                    // FEAST Eigenvalues
  std::vector<ComplexNumber> weights;                            // Quadrature weights.
  const ModelType& model;                                        // DMRG model.
  const LatticeType& lattice;                                    // DMRG lattice.
  ChargeType totalQN;                                            // Overall quantum number associated with the target MPS.
  bool calculateStandardDeviation, screenStandardDeviation;      // If true, calculates the standard deviation for each FEAST state.
  double eMin, eMax;                                             // FEAST integration boundaries.
  double standardDeviationScreening;                             // Screening parameter for the standard deviation
  std::shared_ptr<VectorOfMPSs> currentFEASTMPSs;                // Current approximation to the FEAST eigenvalues;
};

} // namespace FeastHelper

#endif // FEAST_HELPER_CLASS
