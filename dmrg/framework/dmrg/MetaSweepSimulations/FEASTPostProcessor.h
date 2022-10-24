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

  FEASTPostProcessor(int numberOfStates, int numberOfQuadrature, const std::vector<ComplexNumber>& w, const ModelType& inputModel,
                     const LatticeType& inputLattice, BaseParameters& parms)
    : nStates(numberOfStates), nQuad(numberOfQuadrature), weights(w), model(inputModel), lattice(inputLattice), calculateVariance(false)
  {
    energies = std::vector<double>(nStates, 0);
    energiesPrev = std::vector<double>(nStates, 0);
    truncatedEnergy = std::vector<double>(nStates, 0);
    variance = std::vector<double>(nStates, 0);
    totalQN = model.total_quantum_numbers(parms);
    if (parms["feast_calculate_variance"] == "yes")
      calculateVariance = true;
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
    // QR of the overlap
    auto zeroComplex = ComplexNumber(0., 0.);
    ComplexMatrixType leftEigenVectors(nStates, nStates, zeroComplex), rightEigenVectors(nStates, nStates, zeroComplex);
    ComplexVectorType alphaVec(nStates, 0.), betaVec(nStates, zeroComplex);
    alps::numeric::ggev(H, B, alphaVec, betaVec, leftEigenVectors, rightEigenVectors, thresholdForRank_);
    // Retrieves energies
    energiesPrev = energies;
    int rank = 0;
    for (int iState = 0; iState < nStates; iState++) {
      if (std::imag(alphaVec[iState]) > 1.0E-10 || std::imag(betaVec[iState]) > 1.0E-10)
        maquis::cout << " WARNING: Energy of the " << iState << "-th state has a non-negligible imaginary part" << std::endl;
      if (std::abs(betaVec[iState]) > thresholdForRank_) {
        energies[iState] = maquis::real(alphaVec[iState]/betaVec[iState]);
        rank += 1;
      }
      else {
        energies[iState] = 0.;
      }
    }
    // Prints out rank
    maquis::cout << " Number of linearly independent FEAST basis vectors: " << rank << std::endl;
    maquis::cout << std::endl;
    // Final copy of the eigenvectors
    feastEigenVectors = rightEigenVectors;
    // maquis::cout << feastEigenVectors << std::endl;
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
    using VectorOfMPSs = std::vector<MPSType>;
    using MatrixOfMPSs = std::map<std::pair<int, int>, MPSType >;
    MatrixOfMPSs mpsTransformed;
    // Variable definition
    int rank = energies.size();
    auto refNorm = ietl::two_norm(mpsContainer->begin()->second[0]);
    auto result = std::make_shared<VectorOfMPSs>(rank);
    for (auto& iMPS: *mpsContainer)
      iMPS.second[0] /= refNorm;
    //#pragma omp parallel for collapse(2)
    // Actual back-transformation
    for (int iOutput = 0; iOutput < rank; iOutput++) {
      for (int iInput = 0; iInput < nStates; iInput++) {
        for (int iQuad = 0; iQuad < nQuad; iQuad++) {
          MPSType mpsToAdd = mpsContainer->operator[](std::make_pair(iInput, iQuad));
          auto scalingFactor = feastEigenVectors(iInput, iOutput)*weights[iQuad]; // /normVector[iInput];
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
          mpsTransformed[std::make_pair(iOutput, iInput)] = compression::l2r_compress(mpsTransformed[std::make_pair(iOutput, iInput)], mMax, 1.0E-16);
      }
    }
    //#pragma omp parallel for
    for (int iOutput = 0; iOutput < rank; iOutput++) {
      result->operator[](iOutput) = mpsTransformed[std::make_pair(iOutput, 0)];
      for (int iInput = 1; iInput < nStates; iInput++) {
        if (truncEach)
          result->operator[](iOutput) = joinAndTruncate(result->operator[](iOutput), mpsTransformed[std::make_pair(iOutput, iInput)], mMax);
        else
          result->operator[](iOutput) = join(result->operator[](iOutput), mpsTransformed[std::make_pair(iOutput, iInput)]);
        //#pragma omp critical (printEnergy) {
        // std::cout << " Truncated Energy for root " << iOutput << " before truncation = " <<
        //   expval(mpsTransformed[std::make_pair(iOutput, 0)], mpo)/overlap(mpsTransformed[std::make_pair(iOutput, 0)], mpsTransformed[std::make_pair(iOutput, 0)]) << std::endl;
        //}
        //#pragma omp critical (printEnergy)
        //{
        // std::cout << " Truncated Energy for root " << iOutput << " = " <<
        //   expval(mpsTransformed[std::make_pair(iOutput, 0)], mpo)/overlap(mpsTransformed[std::make_pair(iOutput, 0)], mpsTransformed[std::make_pair(iOutput, 0)]) << std::endl;
        //}
        //if (calculateVariance) {
        //    auto norm = overlap(mps_transf(iOutput, 0), mps_transf(iOutput, 0));
        //    auto energy = expval(mps_transf(iOutput, 0), mps_transf(iOutput, 0), mpo)/norm;
        //    auto squaredEnergy = expval_squared(mps_transf(iOutput, 0), mps_transf(iOutput, 0), mpo)/norm;
        //    varianceEnergy[iOutput] = std::sqrt(maquis::real(squaredEnergy) - std::norm(energy));
        //    reducedVariance[iOutput] = std::sqrt(maquis::real(squaredEnergy)) - maquis::real(energy);
        //}
      }
      if (!truncEach)
        result->operator[](iOutput) = compression::l2r_compress(result->operator[](iOutput), mMax, 1.0E-16);
      truncatedEnergy[iOutput] = maquis::real(expval(result->operator[](iOutput), mpo)/norm(result->operator[](iOutput)));
    }
    // If requested, calculates the variance
    if (calculateVariance) {
      for (int iState = 0; iState < variance.size(); iState++)
        variance[iState] = this->getVariance(mpo, result->operator[](iState), mMax);
    }
    return result;
  }

  /** @brief Prints the results of the FEAST calculation */
  void printResults() const {
    maquis::cout << " +----------------------------------------------------------+----------------------+" << std::endl;
    maquis::cout << " |   State    |      Old energy      |      New energy      |   Truncated energy   |" << std::endl;
    maquis::cout << " +----------------------------------------------------------+----------------------+" << std::endl;
    for (int iState = 0; iState < energies.size(); iState++)
        maquis::cout << std::setw(13) << std::internal << iState
                     << std::setw(23) << std::right << std::fixed << std::setprecision(8) << energiesPrev[iState]
                     << std::setw(23) << std::right << std::fixed << std::setprecision(8) << energies[iState]
                     << std::setw(23) << std::right << std::fixed << std::setprecision(8) << truncatedEnergy[iState]
                     << std::endl;
    maquis::cout << " +---------------------------------------------------------------------------------+" << std::endl;
    maquis::cout << std::endl;
    // If requested, prints also the variance
    if (calculateVariance) {
      maquis::cout << " +-----------------------------------+" << std::endl;
      maquis::cout << " |   State    |    Energy variance   |" << std::endl;
      maquis::cout << " +-----------------------------------+" << std::endl;
      for (int iState = 0; iState < energies.size(); iState++) {
        maquis::cout << std::setw(13) << std::internal << iState
                     << std::setw(23) << std::right << std::fixed << std::setprecision(8) << variance[iState]
                     << std::endl;
      }
      maquis::cout << " +-----------------------------------+" << std::endl;
      maquis::cout << std::endl;
    }
  }

  /** @brief Getter for the vibrational energy */
  auto getEnergies() const {
    return energies;
  }

  /** @brief Gets the overall energy variation */
  auto getOverallEnergyVariation() const {
    auto overallSum = std::accumulate(energies.begin(), energies.end(), 0.);
    auto oldOverallSum = std::accumulate(energiesPrev.begin(), energiesPrev.end(), 0.);
    return std::abs(overallSum-oldOverallSum);
  }

private:

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

  /** @brief Variance calculator -- needed to screen the FEAST eigenfunctions */
  auto getVariance(const MPOType& mpo, const MPSType& inputMPS, int mMax) {
    auto traitClass = MPOTimesMPSTraitClass<Matrix, SymmGroup>(inputMPS, model, lattice, totalQN, mMax);
    auto outputMPS = traitClass.applyMPO(mpo);
    auto squaredEnergy = (overlap(outputMPS, outputMPS) + 2.*maquis::real(overlap(outputMPS, inputMPS)*mpo.getCoreEnergy())
                           + overlap(inputMPS, inputMPS)*std::norm(mpo.getCoreEnergy()))/norm(inputMPS);
    auto energy = expval(inputMPS, mpo)/norm(inputMPS);
    return maquis::real(squaredEnergy - std::norm(energy));
  }

  // -- Class members --
  std::shared_ptr<ResultContainerType> mpsContainer;                     // Data structure storing the result of the FEAST linear systems.
  int nStates, nQuad;                                                    // FEAST-specific integer parameters.
  std::vector<double> energies, energiesPrev, truncatedEnergy, variance; // FEAST-specific double parameters.
  ComplexMatrixType feastEigenVectors;                                   // FEAST --> eigenvalues transformation matrix.
  static constexpr int thresholdForRank_ = 1.0E-10;                      // Threshold for rank.
  RealVectorType eigenValues;                                            // FEAST Eigenvalues
  std::vector<ComplexNumber> weights;                                    // Quadrature weights.
  const ModelType& model;                                                // DMRG model.
  const LatticeType& lattice;                                            // DMRG lattice.
  ChargeType totalQN;                                                    // Overall quantum number associated with the target MPS.
  bool calculateVariance;                                                // If true, calculates the variance for each FEAST state.
};

} // namespace FeastHelper

#endif // FEAST_HELPER_CLASS
