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
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"

namespace FeastHelper {

/** @brief Class devoted to the post-processing of the FEAST data */
template <class SymmGroup>
class FEASTPostProcessor {
public:
  using Matrix = cmatrix;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using RealValueType = double;
  using RealMatrixType = alps::numeric::matrix<RealValueType>;
  using RealVectorType = alps::numeric::vector<RealValueType>;
  using ComplexNumber = typename std::complex<RealValueType>;
  using ComplexMatrixType = alps::numeric::matrix<ComplexNumber>;
  using DiagonalMatrixType = typename alps::numeric::associated_real_diagonal_matrix<ComplexMatrixType>::type;
  using ComplexVectorType = alps::numeric::vector<ComplexNumber>;
  using ResultContainerType = std::map<std::pair<int, int>, MPSType>;

  FEASTPostProcessor(int numberOfStates, int numberOfQuadrature, const std::vector<ComplexNumber>& w)
    : nStates(numberOfStates), nQuad(numberOfQuadrature), weights(w)
  {
    energies = std::vector<double>(nStates, 0);
    energiesPrev = std::vector<double>(nStates, 0);
    normVector = std::vector<double>(nStates, 0);
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
    // for (int iThread = 0; iThread < omp_get_max_threads(); iThread++) {
    //     H += Hvec[iThread];
    //     B += Bvec[iThread];
    // }
    // auto overlapDeterminant = calculateDeterminant(B);
    // maquis::cout << std::setprecision(12);
    // maquis::cout << std::scientific;
    // maquis::cout << " Hamiltonian matrix in the FEAST subspace" << std::endl;
    // maquis::cout << H << std::endl;
    // maquis::cout << " Overlap matrix of the FEAST subspace" << std::endl;
    // maquis::cout << B << std::endl;
    // maquis::cout << " Determinant of the overlap matrix " << overlapDeterminant << std::endl;
    // maquis::cout << H(0, 0)/B(0, 0) << std::endl;
    for (int i = 0; i < nStates; i++)
      normVector[i] = std::sqrt(std::real(B(i, i)));
    for (int i = 0; i < nStates; i++) {
      for (int j = 0; j < nStates; j++) {
        H(i, j) /= normVector[i]*normVector[j];
        B(i, j) /= normVector[i]*normVector[j];
      }
    }
    // auto overlapDeterminantAfter = calculateDeterminant(B);
    // std::cout << "Determinant of the overlap matrix after " << overlapDeterminantAfter << std::endl;
    // == Matrix diagonalization ==
    // QR of the overlap
    ComplexMatrixType U, V;
    DiagonalMatrixType S;
    alps::numeric::svd(B, U, V, S);
    int rank = 0;
    for (int iElement = 0; iElement < nStates; iElement++)
      if (std::fabs(S(iElement, iElement)) > thresholdForRank_)
        rank += 1;
    // std::cout << " The FEAST overlap matrix has a rank " << rank << std::endl;
    // cmat_type regularizedInverseSquareRoot = svd.matrixU().block(0, 0, n_states, rank_)*
    //                                          svd.singularValues().head(rank_).array().rsqrt().matrix().asDiagonal();
    ComplexMatrixType regularizedInverseSquareRoot(nStates, rank);
    for (int iRow = 0; iRow < nStates; iRow++)
      for (int iCol = 0; iCol < rank; iCol++)
        regularizedInverseSquareRoot(iRow, iCol) = U(iRow, iCol)/std::sqrt(S(iCol, iCol));
    // cmat_type lowdinHamiltonian = regularizedInverseSquareRoot.adjoint()*H*regularizedInverseSquareRoot;
    ComplexMatrixType tmp(nStates, rank), lowdinHamiltonian(rank, rank);
    gemm(H, regularizedInverseSquareRoot, tmp);
    gemm(adjoint(regularizedInverseSquareRoot), tmp, lowdinHamiltonian);
    // Eigen::SelfAdjointEigenSolver<cmat_type> tmpSolver(lowdinHamiltonian);
    // rvec_type eigenvalues = tmpSolver.eigenvalues().real();
    eigenValues = RealVectorType(rank);
    eigenVectors = ComplexMatrixType(rank, rank);
    // maquis::cout << "Regularized inverse square root" << std::endl;
    // maquis::cout << regularizedInverseSquareRoot << std::endl;
    // maquis::cout << "Lowding Hamiltonian" << std::endl;
    // maquis::cout << lowdinHamiltonian << std::endl;
    lowdinHamiltonian = (lowdinHamiltonian + adjoint(lowdinHamiltonian));
    lowdinHamiltonian /= 2.;
    alps::numeric::heev(lowdinHamiltonian, eigenVectors, eigenValues);
    energiesPrev = energies;
    for (int iState = 0; iState < rank; iState++)
      energies[iState] = eigenValues[iState];
    eigenVectorsRescaled = ComplexMatrixType(rank, rank);
    gemm(regularizedInverseSquareRoot, eigenVectors, eigenVectorsRescaled);
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
    using MatrixOfMPSs = std::map<std::pair<int, int>, MPSType >;
    MatrixOfMPSs mpsTransformed;
    // Variable definition
    int rank = energies.size();
    auto refNorm = ietl::two_norm(mpsContainer->begin()->second[0]);
    std::vector<MPSType> result(rank);
    for (auto& iMPS: *mpsContainer)
      iMPS.second[0] /= refNorm;
    //#pragma omp parallel for collapse(2)
    // Actual back-transformation
    for (int iOutput = 0; iOutput < rank; iOutput++) {
      for (int iInput = 0; iInput < nStates; iInput++) {
        for (int iQuad = 0; iQuad < nQuad; iQuad++) {
          MPSType mpsToAdd = mpsContainer->operator[](std::make_pair(iInput, iQuad));
          auto scalingFactor = eigenVectorsRescaled(iInput, iOutput)*weights[iQuad]/normVector[iInput];
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
      result[iOutput] = mpsTransformed[std::make_pair(iOutput, 0)];
      for (int iInput = 1; iInput < nStates; iInput++) {
        if (truncEach)
          result[iOutput] = joinAndTruncate(result[iOutput], mpsTransformed[std::make_pair(iOutput, iInput)], mMax);
        else
          result[iOutput] = join(result[iOutput], mpsTransformed[std::make_pair(iOutput, iInput)]);
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
        result[iOutput] = compression::l2r_compress(result[iOutput], mMax, 1.0E-16);
    }
    return result;
  }

  /** @brief Prints the results of the FEAST calculation */
  void printResults() const {
    // Select only the ones that lie in user specified interval
    std::cout << " +----------------------------------------------+" << std::endl;
    std::cout << " |   State    |   Old energy   |   New energy   |" << std::endl;
    std::cout << " +----------------------------------------------+" << std::endl;
    for (int iState = 0; iState < energies.size(); iState++)
        std::cout << "  " << std::setw(10) << std::internal << iState << "     "
                  << std::setw(12) << std::right << std::fixed << std::setprecision(3)
                  << energiesPrev[iState] << "     "
                  << energies[iState]     << std::endl;
    std::cout << " +----------------------------------------------+" << std::endl;
    std::cout << std::endl;
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

  // -- Class members --
  std::shared_ptr<ResultContainerType> mpsContainer;
  int nStates, nQuad;
  std::vector<double> energies, energiesPrev, normVector;
  ComplexMatrixType eigenVectors, eigenVectorsRescaled;
  static constexpr int thresholdForRank_ = 1.0E-10;
  RealVectorType eigenValues;
  std::vector<ComplexNumber> weights;
};

} // namespace FeastHelper

#endif // FEAST_HELPER_CLASS
