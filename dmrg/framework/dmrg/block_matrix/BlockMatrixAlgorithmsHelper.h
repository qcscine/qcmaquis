/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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
