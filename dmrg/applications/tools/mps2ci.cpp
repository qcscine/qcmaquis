/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include <cmath>
#include <iterator>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <sys/stat.h>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/sim/matrix_types.h"
#include "dmrg/models/MolecularHamiltonians/util.h"
#include "MPS2CIClass.h"

#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#endif

int main(int argc, char ** argv)
{
  if (argc != 3) {
    maquis::cout << "Usage: mps2ci <mps.h5> <determinants_file> " << std::endl;
    maquis::cout << "See J. Chem. Phys. 126, 244109(2007)" << std::endl;
    exit(1);
  }
  maquis::cout.precision(10);
  std::ifstream param_file(argv[1]);
  if (!param_file) {
    maquis::cerr << "Could not open the mps." << std::endl;
    exit(1);
  }
  // Creates the calculator class
  MPSCICalculator<matrix, grp> calculator(argv[1]);
  // Calculates the overlap
  auto completeness = calculator.calculateOverlap(argv[2]);
  maquis::cout << " Completeness : " << completeness << std::endl;
  return 0;
}

