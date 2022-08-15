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
#include "dmrg/sim/matrix_types.h"

/** @brief Checks that the eigenvalues of a real matrix fulfill the trace property  */
BOOST_AUTO_TEST_CASE(CheckRealHermitianDiagonalization) {
  // Initialization
  int size = 10;
  matrix hamiltonianMatrix(size, size, 0.0);
  alps::numeric::vector<double> eigenValues(size, 0.0);
  // Fills the matrix with the actual data
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      hamiltonianMatrix(i, j) = (i == j) ? static_cast<double>(i) : 0.1;
  alps::numeric::heev(hamiltonianMatrix, eigenValues);
  double expectedTrace = size*(size-1)/2, calculatedTrace = 0.;
  for (int i = 0; i < size; i++)
    calculatedTrace += eigenValues(i);
  BOOST_CHECK_CLOSE(calculatedTrace, expectedTrace, 1.0E-10);
}

/** 
 * @brief Checks by hand that the eigenvalues of a complex Hermitean matrix 
 *        fulfill the eigv. definition (so Mv = lambda*v).
 */
BOOST_AUTO_TEST_CASE(CheckComplexHermitianDiagonalization) {
  // Initialization
  int size = 15;
  std::complex<double> complexZero = std::complex<double>(0., 0.);
  cmatrix hamiltonianMatrix(size, size, complexZero), eigenVectors(size, size, complexZero),
          Mv(size, size, complexZero), lambdaV(size, size, complexZero);
  alps::numeric::vector<double> eigenValues(size, 0.);
  // Fills the matrix with the actual data
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      if (i == j)
        hamiltonianMatrix(i, j) = std::complex<double>(i, 0.);
      else if (i > j)
        hamiltonianMatrix(i, j) = std::complex<double>(2., 3.);
      else
        hamiltonianMatrix(i, j) = std::complex<double>(2., -3.);
  alps::numeric::heev(hamiltonianMatrix, eigenVectors, eigenValues);
  // Checks normalization of eigenvectors
  for (int iEigen = 0; iEigen < size; iEigen++) {
    double norm = 0.;
    for (int iDim = 0; iDim < size; iDim++)
      norm += std::norm(eigenVectors(iEigen, iDim));
    BOOST_CHECK_CLOSE(norm, 1., 1.0E-10);
  }
  // Checks that the definition of eigenvector applies.
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      for (int iJunk = 0; iJunk < size; iJunk++)
        Mv(iRow, iCol) += hamiltonianMatrix(iRow, iJunk)*eigenVectors(iJunk, iCol);
      lambdaV(iRow, iCol) = eigenVectors(iRow, iCol)*eigenValues(iCol);
    }
  }
  // Checks matrix equality element-wise
  for (int iRow = 0; iRow < size; iRow++)
    for (int iCol = 0; iCol < size; iCol++)
      BOOST_CHECK_CLOSE(std::norm(lambdaV(iRow, iCol)), std::norm(Mv(iRow, iCol)), 1.0E-10);
}

/** 
 * @brief Checks by hand that the eigenvalues of a complex non Hermitean matrix 
 *        fulfill the eigv. definition (so Mv = lambda*v).
 */
BOOST_AUTO_TEST_CASE(CheckComplexNonHermitianDiagonalization) {
  // Initialization
  int size = 20;
  std::complex<double> complexZero = std::complex<double>(0., 0.);
  cmatrix hamiltonianMatrix(size, size, complexZero), eigenVectorsLeft(size, size, complexZero),
          eigenVectorsRight(size, size, complexZero), Mv(size, size, complexZero),
          lambdaV(size, size, complexZero);
  alps::numeric::vector<std::complex<double>> eigenValues(size, complexZero);
  // Fills the matrix with the actual data
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      if (i == j)
        hamiltonianMatrix(i, j) = std::complex<double>(i, size-i);
      else if (i > j)
        hamiltonianMatrix(i, j) = std::complex<double>(4., 3.);
      else
        hamiltonianMatrix(i, j) = std::complex<double>(2., 3.);
  alps::numeric::geev(hamiltonianMatrix, eigenVectorsLeft, eigenVectorsRight, eigenValues);
  // Does the check for the right eigenvectors.
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      for (int iJunk = 0; iJunk < size; iJunk++)
        Mv(iRow, iCol) += hamiltonianMatrix(iRow, iJunk)*eigenVectorsRight(iJunk, iCol);
      lambdaV(iRow, iCol) = eigenVectorsRight(iRow, iCol)*eigenValues(iCol);
    }
  }
  // Checks matrix equality element-wise
  for (int iRow = 0; iRow < size; iRow++)
    for (int iCol = 0; iCol < size; iCol++)
      BOOST_CHECK_CLOSE(std::norm(lambdaV(iRow, iCol)), std::norm(Mv(iRow, iCol)), 1.0E-10);
}