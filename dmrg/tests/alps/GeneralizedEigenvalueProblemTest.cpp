/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE alps

#include <random>
#include <boost/mpl/assert.hpp>
#include <boost/test/included/unit_test.hpp>
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

/** @brief Checks that the left and right eigenvectors are correct */
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

/** @brief Checks that the generalized eigenvalue solver works for complex matrices */
BOOST_AUTO_TEST_CASE(CheckGeneralizedEigenvalueComplex) {
  int size = 10;
  cmatrix complexHamiltonian(size, size), leftEigenVectors(size, size, std::complex<double>(0., 0.)),
    rightEigenVectors(size, size, std::complex<double>(0., 0.)), complexOverlap(size, size, std::complex<double>(0., 0.));
  alps::numeric::vector<std::complex<double>> eigenValues(size, std::complex<double>(0., 0.));
  // Constructs the random engine.
  std::mt19937 gen(1991);
  std::uniform_real_distribution<double> diagonalDistribution(0.1, 1.), offDiagonalDistribution(0., 0.1);
  for (int iRow = 0; iRow < size; iRow++) {
    complexHamiltonian(iRow, iRow) = std::complex<double>(diagonalDistribution(gen), 0.);
    complexOverlap(iRow, iRow) = std::complex<double>(diagonalDistribution(gen), 0.);
    for (int iCol = 0; iCol <= iRow; iCol++) {
      complexHamiltonian(iRow, iCol) = std::complex<double>(diagonalDistribution(gen), diagonalDistribution(gen));
      complexOverlap(iRow, iCol) = std::complex<double>(diagonalDistribution(gen), diagonalDistribution(gen));
      complexHamiltonian(iCol, iRow) = std::conj(complexHamiltonian(iRow, iCol));
      complexOverlap(iCol, iRow) = std::conj(complexOverlap(iRow, iCol));
    }
  }
  alps::numeric::ggev(complexHamiltonian, complexOverlap, eigenValues, leftEigenVectors, rightEigenVectors);
  // Verifies that the right eigenvalues are correct
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      std::complex<double> tmp = 0.;
      for (int iDummy = 0; iDummy < size; iDummy++) {
        tmp += complexHamiltonian(iRow, iDummy)*rightEigenVectors(iDummy, iCol);
        tmp -= complexOverlap(iRow, iDummy)*rightEigenVectors(iDummy, iCol)*eigenValues[iCol];
      }
      BOOST_CHECK_SMALL(std::abs(tmp), 1.0E-10);
    }
  }
  // Verifies that the left eigenvalues are correct
  for (int iRow = 0; iRow < size; iRow++) {
    for (int iCol = 0; iCol < size; iCol++) {
      std::complex<double> tmp = 0.;
      for (int iDummy = 0; iDummy < size; iDummy++) {
        tmp += std::conj(leftEigenVectors(iDummy, iRow))*complexHamiltonian(iDummy, iCol);
        tmp -= eigenValues[iRow]*std::conj(leftEigenVectors(iDummy, iRow))*complexOverlap(iDummy, iCol);
      }
      BOOST_CHECK_SMALL(std::abs(tmp), 1.0E-10);
    }
  }
}