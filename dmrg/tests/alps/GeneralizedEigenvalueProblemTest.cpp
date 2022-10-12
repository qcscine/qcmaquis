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

#define BOOST_TEST_MODULE alps

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/block_matrix/detail/alps.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "dmrg/sim/matrix_types.h"

/** @brief Checks consistency between heef and ggev for a conventional eigenvalue problem */
BOOST_AUTO_TEST_CASE(CheckGeneralizedEigenvalueProblemTrivial) {
  int size = 10;
  matrix hamiltonianMatrix(size, size, 0.0), leftEigenVectors(size, size, 0.0), 
    rightEigenVectors(size, size, 0.), overlap(size, size, 0.0);
  alps::numeric::vector<std::complex<double> > eigenValues(size, 0.0);
  alps::numeric::vector<double> eigenValuesReal(size, 0.0);
  // Fills the matrix with the actual data
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      hamiltonianMatrix(i, j) = (i == j) ? static_cast<double>(i) : 0.1;
    overlap(i, i) = 1.;
  }
  auto hamiltonianMatrixCopy = hamiltonianMatrix;
  // GGEV
  alps::numeric::ggev(hamiltonianMatrix, overlap, eigenValues, leftEigenVectors, rightEigenVectors);
  std::sort(eigenValues.begin(), eigenValues.end(), [](auto a, auto b) { return std::real(a) > std::real(b); });
  // HHEV
  alps::numeric::heev(hamiltonianMatrixCopy, eigenValuesReal);
  for (int iElement = 0; iElement < size; iElement++)
    BOOST_CHECK_CLOSE(std::real(eigenValues[iElement]), eigenValuesReal[iElement], 1.0E-10);
}

BOOST_AUTO_TEST_CASE(CheckGeneralizedEigenvalueLeftAndRight) {
  int size = 2;
  matrix hamiltonianMatrix(size, size, 0.0), leftEigenVectors(size, size, 0.0), 
    rightEigenVectors(size, size, 0.), overlap(size, size, 0.0);
  alps::numeric::vector<std::complex<double> > eigenValues(size, 0.0);
  hamiltonianMatrix(0, 0) = -2.0;
  hamiltonianMatrix(1, 1) = -1.0;
  hamiltonianMatrix(1, 0) =  0.2;
  hamiltonianMatrix(0, 1) =  0.3;
  overlap(0, 0) = 1.0;
  overlap(1, 1) = 1.0;
  overlap(0, 1) = 0.2;
  overlap(1, 0) = 0.2;
  alps::numeric::ggev(hamiltonianMatrix, overlap, eigenValues, leftEigenVectors, rightEigenVectors);
  // Verifies that the right eigenvalues are correct
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      double tmp = 0.;
      for (int iDummy = 0; iDummy < size; iDummy++) {
        tmp += hamiltonianMatrix(iRow, iDummy)*rightEigenVectors(iDummy, iCol);
        tmp -= overlap(iRow, iDummy)*rightEigenVectors(iDummy, iCol)*std::real(eigenValues[iCol]);
      }
      BOOST_CHECK_SMALL(tmp, 1.0E-10);
    }
  }
  // Verifies that the left eigenvalues are correct
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      double tmp = 0.;
      for (int iDummy = 0; iDummy < size; iDummy++) {
        tmp += leftEigenVectors(iDummy, iRow)*hamiltonianMatrix(iDummy, iCol);
        tmp -= std::real(eigenValues[iRow])*leftEigenVectors(iDummy, iRow)*overlap(iDummy, iCol);
      }
      BOOST_CHECK_SMALL(tmp, 1.0E-10);
    }
  }
}