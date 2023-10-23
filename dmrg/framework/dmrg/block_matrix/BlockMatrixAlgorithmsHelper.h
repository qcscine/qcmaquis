/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef BLOCK_MATRIX_ALGORITHMS_HELPER_H
#define BLOCK_MATRIX_ALGORITHMS_HELPER_H

#include "dmrg/sim/matrix_types.h"

/** @brief Helper class for block matrix algorithms */
template<class Matrix, class SymmGroup>
class BlockMatrixAlgorithmsHelperClass {
public:
  /**
   * @brief Adjusts the phase of the
   */
  static void adjustPhase(Matrix& eigenVectors) {};
};

/** @brief Specialization for complex-values matrices */
template<class SymmGroup>
class BlockMatrixAlgorithmsHelperClass<cmatrix, SymmGroup> {
  using ComplexType = std::complex<double>;
public:
  static void adjustPhase(cmatrix& eigenVectors) {
    std::vector<ComplexType> vectorOfPhases;
    // Picks up the phases
    for (int iCol = 0; iCol < num_cols(eigenVectors); iCol++) {
      int iRow = 0;
      bool exit = false;
      do {
        if (std::abs(eigenVectors(iRow, iCol)) > threshold) {
          vectorOfPhases.push_back(std::exp(-std::complex<double>(0., 1.)*std::arg(eigenVectors(iRow, iCol))));
          exit = true;
        }
        else if (iRow == num_rows(eigenVectors)-1) {
          vectorOfPhases.push_back(std::complex<double>(1., 0.));
          exit = true;
        }
        else {
          iRow++;
        }
      } while (!exit);
    }
    // Final rescaling
    for (int iCol = 0; iCol < num_cols(eigenVectors); iCol++)
      for (int iRow = 0; iRow < num_rows(eigenVectors); iRow++)
        eigenVectors(iRow, iCol) *= vectorOfPhases[iCol];
  }
private:
  static constexpr double threshold = 1.0E-16;
};

#endif // BLOCK_MATRIX_ALGORITHMS_HELPER_H
