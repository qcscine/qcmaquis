/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE mpstensor

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>

#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/BenzeneFixture.h"

/** @brief Test for the constructor from an eigen tensor */
BOOST_FIXTURE_TEST_CASE(MPSTensorFromEigen, BenzeneFixture) {
#ifdef HAVE_TwoU1PG
  using SymmetryType = TwoU1PG;
  using ModelType = Model<matrix, SymmetryType>;
  using MPSTensorType = MPSTensor<matrix, SymmetryType>;
  // Generates the HF MPS
  auto lattice = Lattice(parametersBenzene);
  auto modelHF = Model<matrix, SymmetryType>(lattice, parametersBenzene);
  auto mpsHF = MPS<matrix, SymmetryType>(lattice.size(), *(modelHF.initializer(lattice, parametersBenzene)));
  // Just picks randomly the third site
  auto mpsTensor = mpsHF[3];
  auto normFromMPS = mpsTensor.scalar_norm();
  auto normFromEigen = std::sqrt(mpsTensor.getEigenRepresentation().squaredNorm());
  BOOST_CHECK_CLOSE(normFromMPS, normFromEigen, 1.0E-10);
#endif // HAVE_TwoU1PG
}

#ifdef HAVE_TwoU1PG
BOOST_FIXTURE_TEST_CASE(EigenFromMPSTensor, BenzeneFixture) {
  using SymmetryType = TwoU1PG;
  using ModelType = Model<matrix, SymmetryType>;
  using MPSTensorType = MPSTensor<matrix, SymmetryType>;
  // Generates the HF MPS
  auto lattice = Lattice(parametersBenzene);
  auto modelHF = Model<matrix, SymmetryType>(lattice, parametersBenzene);
  auto mpsHF = MPS<matrix, SymmetryType>(lattice.size(), *(modelHF.initializer(lattice, parametersBenzene)));
  // Picks now the 2-nd mps tensor
  auto mpsTensor = mpsHF[2];
  auto numElements = mpsTensor.num_elements();
  Eigen::Matrix<double, Eigen::Dynamic, 1> vectorInitialize = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(numElements);
  vectorInitialize[0] = 3.;
  vectorInitialize[5] = 4.;
  mpsTensor.fillWithEigenVector(vectorInitialize);
  auto finalNorm = mpsTensor.scalar_norm();
  BOOST_CHECK_CLOSE(finalNorm, 5., 1.0E-14);
}
#endif // HAVE_TwoU1PG