/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#include <iostream>
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/compression.h"
#include "dmrg/models/MolecularHamiltonians/transform_symmetry.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>

#include "copydir.hpp"

#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#elif defined(USE_SU2U1)
typedef SU2U1 grp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#endif

template <class Matrix, class SymmGroup>
void compress_mps(MPS<Matrix, SymmGroup>& mps, std::size_t Mmax) {
  maquis::cout << "Compressing MPS, new m_max = " << Mmax << std::endl;

  matrix::value_type final_norm = norm(mps);
  matrix::value_type compression_trace = 1.0;

  mps = compression::l2r_compress(mps, Mmax, 1e-8, compression_trace);
  maquis::cout << "- compression trace          : " << compression_trace
               << std::endl;
  mps[0].multiply_by_scalar(compression_trace * sqrt(final_norm));
}

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      std::cout << "Usage: " << argv[0] << " mps_checkpoint.h5 new_m"
                << std::endl;
      return 1;
    }

    std::size_t m;
    std::string argv1 = argv[1];
    // read in the m value
    try {
      m = boost::lexical_cast<std::size_t>(argv[2]);
    } catch (std::exception& e) {
      std::cerr << "Cannot read the m value!" << std::endl;
      return 1;
    }

    MPS<matrix, grp> mps;
    load(argv1, mps);
    // backup the original mps
    std::string backupName =
        boost::ireplace_first_copy(argv1, ".h5", "_uncompressed.h5");
    copyDir(argv1, backupName);
    maquis::cout << "Original MPS saved on " << backupName << std::endl;

    compress_mps(mps, m);
    save(argv1, mps);

  } catch (std::exception& e) {
    std::cerr << "Error:" << std::endl << e.what() << std::endl;
    return 1;
  }
}
